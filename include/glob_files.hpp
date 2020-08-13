/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef GLOB_H_GUARD
#define GLOB_H_GUARD
#include <glob.h>
#include <string.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::string> glob(const std::string& pattern);

#endif
