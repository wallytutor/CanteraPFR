// ***************************************************************************
// Provides utilities for parsing XML.
//
// Author : Walter Dal'Maz Silva
// Date   : December 24th 2018
// ***************************************************************************

#ifndef __XML_XMLPARSE_HPP__
#define __XML_XMLPARSE_HPP__

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

// Auxiliary error/warning messages.
#define DEBUG_FMT "\n[DEBUG] <%s(%s)> l.%d\n"
#define ERROR_FMT "\n[ERROR] <%s(%s)> : "
#define WARNS_FMT "\n[WARNING] <%s(%s)> : "
#define MACRO_DEBUG fprintf(stderr, DEBUG_FMT,__FILE__,__func__,__LINE__)
#define MACRO_ERROR fprintf(stderr, ERROR_FMT,__FILE__,__func__)
#define MACRO_WARNS fprintf(stderr, WARNS_FMT,__FILE__,__func__)

namespace bpt = boost::property_tree;

namespace Cantera
{

constexpr char TAB[] = "    ";

/*! \brief Print data formatted in JSON-like format.
 *
 *  Provides recursive printing of a data tree consisting of many levels. Each
 *  penetration into a deeper level is done by providing an extra tabulation.
 *
 *  \param pt property tree to produce output.
 *  \param tabulation tabulation level of output.
 */
inline void printJSON(const bpt::ptree& pt, unsigned tabulation=0)
{
    for (unsigned i = 0; i != tabulation; ++i) { std::cout << Cantera::TAB; }

    for (auto it = pt.begin(); it != pt.end(); ++it)
    {
        std::string st = it->first;
        std::string nd = it->second.get_value<std::string>();
        std::cout << st << " : " << nd << std::endl;
        printJSON(it->second, tabulation+1);
    }
}

/*! \brief Retrieves arguments from property tree node or attribute.
 *
 *  Access to argument in a property tree is provided by key access. The method
 *  implements formatted error output handling and a default argument can be
 *  returned when presence in tree is not mandatory.
 *
 *  \param pt property tree to retrieve node or attribute value.
 *  \param key key to node or attribute to retrive, including multi-level.
 *  \param msg error message to be print or raised if necessary.
 *  \param defval default value of same type as return object.
 *  \param required whether or not the key has to be present.
 */
template<class T>
T get_argument(const bpt::ptree& pt, const std::string& key,
    const std::string& msg, T defval, bool required=false)
{
    try
    {
        return pt.get<T>(key);
    }
    catch (std::exception& err)
    {
        if (!required)
        {
            return defval;
        }
        MACRO_ERROR;
        throw std::runtime_error(err.what());
    }
}

/*! \brief Reads a whole file into a string stream.
 *
 *  File path passed as argument to this function is read into a string buffer
 *  for later processing. No existance test is performed and an error is thrown
 *  in case it was not possible to read the file.
 *
 *  \param fpath path to file to read into string buffer.
 */
inline std::stringstream readFile(const std::string& fpath)
{
    try
    {
        std::ifstream ssfile(fpath);

        if (ssfile)
        {
            std::stringstream buffer;
            buffer << ssfile.rdbuf();
            ssfile.close();
            return buffer;
        }

        throw std::runtime_error("Unable to read " + fpath);
    }
    catch (std::exception& err)
    {
        MACRO_ERROR;
        throw std::runtime_error(err.what());
    }
}

/*! Reads a XML file into a property tree.
 *
 *  A property tree is loaded from a string buffer. No file of path checks are
 *  performed and an error is raised if loading is not possible.
 *
 *  \param pt property tree to load data into.
 *  \param fpath path to file to read XML data into property tree.
 *  \param option XML parsing option, default 2 strips comments.
 */
inline void readXMLFile(bpt::ptree& pt, const std::string& fpath, int opt = 2)
{
    try
    {
        std::stringstream buffer = readFile(fpath);
        bpt::read_xml(buffer, pt, opt);
    }
    catch (std::exception& err)
    {
        MACRO_ERROR;
        throw std::runtime_error(err.what());
    }
}

} // (namespace Cantera)

#endif // (__XML_XMLPARSE_HPP__)

// ###########################################################################
//         EOF
// ###########################################################################
