/**
 * @file VectorFile.hpp
 * @brief Base abstract class for vector files.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTORFILE_HPP
#define WILDRIVER_VECTORFILE_HPP




#include <string>
#include <memory>

#include "IVectorFile.hpp"
#include "Vector.hpp"




namespace WildRiver
{


class VectorFile :
  public IVectorFile

{
  public:
    /**
     * @brief Default constructor which intializes the vector properties. 
     */
    VectorFile();


    /**
     * @brief Virtual destructor.
     */
    virtual ~VectorFile();


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    virtual ind_t getSize() override
    {
      return m_vector.getSize();
    }


    /**
     * @brief Set the size of the vector.  
     *
     * @param size The new size of the vector.
     */
    virtual void setSize(
        const ind_t size) override
    {
      m_vector.setSize(size);
    }


    /**
     * @brief Check to see if the size of the vector has been set.
     *
     * @return True if the size has been set, false otherwise. 
     */
    bool isSizeSet() const noexcept
    {
      return m_vector.isSizeSet();
    }


  private:
    /**
     * @brief Underlying vector structure.
     */
    Vector m_vector;




};




}




#endif
