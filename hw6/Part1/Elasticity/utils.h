/* CS/CNS 171
 * Fall 2015
 * Written by Kevin (Kevli) Li (Class of 2016)
 *
 * Straightforward utility functions.
 */

#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <string>
#include <vector>

/* Function prototypes */

/* Absolute value function */
static float abs(float n);

/* Splits given string based on given delimiter and stores tokens in given vector */
static std::vector<std::string> &split(const std::string &str,
                                       char delim,
                                       std::vector<std::string> &tokens);
/* Splits given string based on given delimiter and returns tokens in a vector */
static std::vector<std::string> split(const std::string &str, char delim);

/* Function implementations */

static float abs(float n)
{
    return (n < 0) ? n * -1.0 : n;
}

static std::vector<std::string> &split(const std::string &str,
                                       char delim,
                                       std::vector<std::string> &tokens)
{
    std::stringstream ss(str);
    std::string token;

    while (getline(ss, token, delim))
        if(!token.empty())
            tokens.push_back(token);

    return tokens;
}

static std::vector<std::string> split(const std::string &str, char delim)
{
    std::vector<std::string> tokens;
    split(str, delim, tokens);
    return tokens;
}

#endif
