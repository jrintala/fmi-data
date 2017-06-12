#include "generator.h"

namespace himan{

/*
 * Finds the maximum values for a series of info_t on the interval [begin,end)
 * */

template <class InputIt>
std::shared_ptr<himan::info> Max(InputIt begin, InputIt end)
{
	// If time series has zero length
	if (begin == end) return nullptr;

	// Find first field that cntains data
        while (*begin == nullptr)
        {
                ++begin;
                if (begin == end)
                        return nullptr;
        }

	// Set first field as first set of maximum values
	auto maxInfo = *begin;
	maxInfo->ReGrid();
	++begin;

	// Update set of maximum values
	for (; begin != end; ++begin)
	{
                if (*begin ==  nullptr)
                        continue;

		auto in = VEC((*begin)).begin();
		auto out = VEC(maxInfo).begin();
		for (; in != VEC((*begin)).end(), out != VEC(maxInfo).end(); ++in, ++out)
		{
			*out = std::max(*in, *out);
		}
	}

	return maxInfo;
}

/*
 * Finds the minimum values for a series of info_t on the interval [begin,end)
 * */

template <class InputIt>
std::shared_ptr<himan::info> Min(InputIt begin, InputIt end)
{
        if (begin == end) return nullptr;

	while (*begin == nullptr)
	{
		++begin;
		if (begin == end)
			return nullptr;
	}

        auto minInfo = *begin;
        minInfo->ReGrid();

        for (; begin != end; ++begin)
        {
		if (*begin ==  nullptr)
			continue;

                auto in = VEC((*begin)).begin();
                auto out = VEC(minInfo).begin();
                for (; out != VEC((*begin)).end(), out != VEC(minInfo).end(); ++in, ++out)
                {
                        *out = std::min(*in, *out);
                }
        }

        return minInfo;
}

}
