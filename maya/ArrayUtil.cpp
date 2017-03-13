#include "ArrayUtil.h"

namespace ArrayUtil
{
	int* rotateArrayCopy(int* src, int size, int amount)
	{
		int* dst = new int[size];
		for (int sidx = amount, didx = 0; sidx < size; sidx++, didx++)
		{
			dst[didx] = src[sidx];
		}
		for (int sidx = 0, didx = size - amount; didx < size; sidx++, didx++)
		{
			dst[didx] = src[sidx];
		}
		return dst;
	}

	int* range(int rangeStart, int rangeEnd)
	{
		int size = rangeEnd - rangeStart;
		int* arr = new int[size];
		for (int i = 0; i < size; i++)
		{
			arr[i] = rangeStart + i;
		}
		return arr;
	}
}
