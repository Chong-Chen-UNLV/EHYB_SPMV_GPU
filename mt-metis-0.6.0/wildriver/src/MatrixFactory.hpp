/**
 * @file MatrixFactory.hpp
 * @brief Class for instantiating matrix files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-05
 */




#ifndef WILDRIVER_MATRIXFACTORY_HPP
#define WILDRIVER_MATRIXFACTORY_HPP




#include <memory>
#include <string>
#include "IMatrixFile.hpp"



namespace WildRiver
{


class MatrixFactory
{
  public:
    /**
     * @brief Instantiate a matrix file.
     *
     * @param fname The name of the file to open.
     *
     * @return The instantiated object.
     *
     * @throw UnknownExtensionException If no class reading the specified file
     * type can be found.
     */
    static std::unique_ptr<IMatrixFile> make(
        std::string const & name);
};




}




#endif
