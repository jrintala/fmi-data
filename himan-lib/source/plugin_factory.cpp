/*
 * plugin_factory.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: partio
 */

#include "plugin_factory.h"
#include "logger_factory.h"
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <cstdlib>
#include <dlfcn.h>
#include "util.h"

using namespace himan;
using namespace himan::plugin;

plugin_factory* plugin_factory::itsInstance = NULL;

plugin_factory* plugin_factory::Instance()
{
    if (!itsInstance)
    {
        itsInstance = new plugin_factory();
    }

    return itsInstance;
}


plugin_factory::plugin_factory() : itsPluginSearchPath()
{

    itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("plugin_factory"));

    char* path;

    path = std::getenv("HIMAN_LIBRARY_PATH");

    if (path != NULL)
    {
        std::vector<std::string> paths = util::Split(std::string(path), ":", false);

        itsPluginSearchPath.insert(itsPluginSearchPath.end(), paths.begin(), paths.end());

    }
    else
    {
        itsLogger->Trace("Environment variable HIMAN_LIBRARY_PATH not set -- search plugins only from pre-defined locations");
    }

    itsPluginSearchPath.push_back(".");
    itsPluginSearchPath.push_back("/usr/lib64/himan-plugins"); // Default RPM location

    ReadPlugins();
}

// Hide constructor
plugin_factory::~plugin_factory() {}

std::vector<std::shared_ptr<himan_plugin> > plugin_factory::Plugins(HPPluginClass pluginClass)
{

    std::vector<std::shared_ptr<himan_plugin>> thePlugins;

    for (size_t i = 0; i < itsPluginFactory.size(); i++)
    {
        if (pluginClass == itsPluginFactory[i]->Plugin()->PluginClass())
        {
            thePlugins.push_back(itsPluginFactory[i]->Plugin());
        }

        else if (pluginClass == kUnknownPlugin)
        {
            thePlugins.push_back(itsPluginFactory[i]->Plugin());
        }
    }

    return thePlugins;
}

std::vector<std::shared_ptr<himan_plugin> > plugin_factory::CompiledPlugins()
{
    return Plugins(kCompiled);
}
std::vector<std::shared_ptr<himan_plugin> > plugin_factory::AuxiliaryPlugins()
{
    return Plugins(kAuxiliary);
}
std::vector<std::shared_ptr<himan_plugin> > plugin_factory::InterpretedPlugins()
{
    return Plugins(kInterpreted);
}

/*
 * Plugin()
 *
 * Return instance of the requested plugin if found. Caller must cast
 * the plugin to the derived class. If second argument is true, a new
 * instance is created and returned. Otherwise function behaves like
 * a regular factory pattern and return one known instance to each
 * caller (this is suitable only in non-threaded functions).
 */

std::shared_ptr<himan_plugin> plugin_factory::Plugin(const std::string& theClassName, bool theNewInstance)
{

    for (size_t i = 0; i < itsPluginFactory.size(); i++)
    {

        if ((itsPluginFactory[i]->Plugin()->ClassName() == theClassName) ||
                (itsPluginFactory[i]->Plugin()->ClassName() == "himan::plugin::" + theClassName))
        {
            if (theNewInstance)
            {
                return itsPluginFactory[i]->Clone();
            }
            else
            {
                return itsPluginFactory[i]->Plugin();
            }
        }
    }

    throw std::runtime_error("plugin_factory: Unknown plugin clone operation requested: " + theClassName);
}

/*
 * ReadPlugins()
 *
 * Read plugins from defined paths. Will try to load all files in given directories
 * that end with .so. Will not ascend to child directories (equals to "--max-depth 1").
 */

void plugin_factory::ReadPlugins()
{

    using namespace boost::filesystem;

    directory_iterator end_iter;

    for (size_t i = 0; i < itsPluginSearchPath.size(); i++)
    {
    	itsLogger->Trace("Search plugins from " + itsPluginSearchPath[i]);

        path p (itsPluginSearchPath[i]);

        try
        {
            if (exists(p) && is_directory(p))      // is p a directory?
            {

                for ( directory_iterator dir_iter(p) ; dir_iter != end_iter ; ++dir_iter)
                {
                    if (dir_iter->path().filename().extension().string() == ".so")
                    {
                        Load(dir_iter->path().string());
                    }
                }
            }
        }

        catch (const filesystem_error& ex)
        {
            itsLogger->Error(ex.what());
        }
    }

}

bool plugin_factory::Load(const std::string& thePluginFileName)
{

    /*
     * Open libraries with
     *
     *   RTLD_LAZY
     * We don't specify the ordering of the plugin load -- usually it is alphabetical
     * but not necessarily so. With RTLD_LAZY the symbols aren't checked during load
     * which means that the first loaded plugin can refer to functions defined in the
     * last loaded plugin without compiler issuing warnings like
     *  <file>.so undefined symbol: ...
     *
     *   RTLD_GLOBAL
     * We need this because core library (himan-lib) need to access aux plugins
     * and the plugin symbol information needs to be propagated throughout the
     * plugins system (or using aux plugins in core lib will fail).
     */

    void* theLibraryHandle = dlopen(thePluginFileName.c_str(), RTLD_LAZY | RTLD_GLOBAL);

    if (!theLibraryHandle)
    {
        itsLogger->Error("Unable to load plugin: " + std::string(dlerror()));
        return false;
    }

    dlerror(); // clear error handle

    //create_t* create_plugin = (create_t*) dlsym(theLibraryHandle, "create");
    create_t* create_plugin =  reinterpret_cast<create_t*> (dlsym(theLibraryHandle, "create"));

    if (!create_plugin)
    {
        itsLogger->Error("Unable to load symbol: " + std::string(dlerror()));
        return false;
    }

    std::shared_ptr<plugin_container> mc = std::shared_ptr<plugin_container> (new plugin_container(theLibraryHandle, create_plugin()));

    for (size_t i = 0; i < itsPluginFactory.size(); i++)
    {
        if (mc->Plugin()->ClassName() == itsPluginFactory[i]->Plugin()->ClassName())
        {
            itsLogger->Trace("Plugin '" + mc->Plugin()->ClassName() + "' found more than once, skipping one found from '" + thePluginFileName + "'");
            return true;
        }
    }

    itsPluginFactory.push_back(mc);

    itsLogger->Debug("Load " + thePluginFileName);

    return true;
}
