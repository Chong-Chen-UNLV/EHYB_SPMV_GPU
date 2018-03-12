/**
 * @file Vector.hpp
 * @brief Class for vectors.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 * @date 2016-02-07
 */




#ifndef WILDRIVER_VECTOR_HPP
#define WILDRIVER_VECTOR_HPP



namespace WildRiver
{


class Vector
{
  public:
    /**
     * @brief Constant value for an unset size.
     */
    static constexpr ind_t UNSET_SIZE = static_cast<ind_t>(-1);


    /**
     * @brief Create a new vector with an unset size.
     */
    Vector() noexcept :
      m_size(UNSET_SIZE)
    {
      // do nothing
    }


    /**
     * @brief Virtual destructor.
     */
    ~Vector()
    {
      // do nothing
    }


    /**
     * @brief Get the size of the vector. May alter the internal state of the
     * reader.
     *
     * @return The size of the vector.
     */
    ind_t getSize() const noexcept
    {
      return m_size;
    }


    /**
     * @brief Set the size of the vector.  
     *
     * @param size The new size of the vector.
     */
    void setSize(
        const ind_t size) noexcept
    {
      m_size = size;
    }


    /**
     * @brief Check to see if the size of the vector has been set.
     *
     * @return True if the size has been set, false otherwise. 
     */
    bool isSizeSet() const noexcept
    {
      return m_size != UNSET_SIZE; 
    }


  private:
    /**
     * @brief The size of the vector.
     */
    ind_t m_size;




};




}




#endif
