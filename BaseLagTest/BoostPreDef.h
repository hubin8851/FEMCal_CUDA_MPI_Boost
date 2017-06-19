#pragma once
//基本头文件
#define  BOOST_ALL_NO_LIB
#include <boost/enable_shared_from_this.hpp>
//#include <boost/timer.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/system/system_error.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

//字符串相关操作头文件
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
//文件操作头文件
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
//定时器头文件
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
//线程池头文件
#include <boost/thread.hpp>
#include <boost/threadpool/pool.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread_pool.hpp>
//XML方面的
#include <boost/property_tree/ptree.hpp>  
#include <boost/property_tree/xml_parser.hpp>  
#include <boost/typeof/typeof.hpp>
//cpu_timer高精度计时器
#include <boost/timer/timer.hpp>

#include <boost/bimap.hpp>