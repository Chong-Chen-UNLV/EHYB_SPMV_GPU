/**
 * @file VectorFactory.hpp
 * @brief Class for instantiating vector files. 
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTORFACTORY_HPP
#define WILDRIVER_VECTORFACTORY_HPP




#include <memory>
#include <string>

#include "IVectorFile.hpp"




namespace WildRiver
{


class VectorFactory
{
  public:
    /**
     * @brief Allocate a new vector file subclass based on teh file extension.
     *
     * @param name The filename/path to open.
     *
     * @return The newly opened vector file.
     */
    static std::unique_ptr<IVectorFile> make(
        std::string const & name);


};




}




#endif
