#pragma once
#include <string>
#include <iostream>
#include <fstream>

class csvfile;

inline static csvfile& endrow(csvfile& file);
inline static csvfile& flush(csvfile& file);

class csvfile
{
	std::ofstream fs_;
	const std::string separator_;
    std::string filename_string;
public:
    csvfile() { }
	csvfile(const std::string filename, const std::string separator = ";")
		: fs_()
		, separator_(separator)
	{
		fs_.exceptions(std::ios::failbit | std::ios::badbit);
		fs_.open(filename);
        filename_string = filename;
	}

	~csvfile()
	{
		flush();
		fs_.close();
	}

    std::string getFileName()
    {
        return filename_string;
    }

	void flush()
	{
		fs_.flush();
	}

	void endrow()
	{
		fs_ << std::endl;
	}

    void precision(int n)
    {
        fs_.precision(n);
    }

	csvfile& operator << ( csvfile& (* val)(csvfile&))
	{
		return val(*this);
	}

	csvfile& operator << (const char * val)
	{
		fs_ << '"' << val << '"' << separator_;
		return *this;
	}

	csvfile& operator << (const std::string & val)
	{
		fs_ << '"' << val << '"' << separator_;
		return *this;
	}

	template<typename T>
	csvfile& operator << (const T& val)
	{
		fs_ << val << separator_;
		return *this;
	}
};


inline static csvfile& endrow(csvfile& file)
{
	file.endrow();
	return file;
}

inline static csvfile& flush(csvfile& file)
{
	file.flush();
	return file;
}

