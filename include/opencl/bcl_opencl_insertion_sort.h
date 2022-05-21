// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_OPENCL_INSERTION_SORT_H_
#define BCL_OPENCL_INSERTION_SORT_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "linal/bcl_linal_matrix.h"
#include "model/bcl_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InsertionSort
    //! @brief performs insertion sort on columns of matrix for k elements of interest also providing index matrix
    //!
    //! @see @link example_opencl_insertion_sort.cpp @endlink
    //! @author loweew
    //! @date Mar 22, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API InsertionSort :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! index matrix on device
      mutable Matrix< int>      m_IndexMatrixOnDevice;

      //! opencl queue
      CommandQueue              m_Queue;

      //! opencl program
      cl::Program               m_Program;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from command queue
      InsertionSort();

      //! @brief constructor from command queue
      //! @param QUEUE command queue
      InsertionSort( const CommandQueue &QUEUE);

      //! @brief Clone function
      //! @return pointer to new InsertionSort
      InsertionSort *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get index matrix
      //! @return index matrix
      Matrix< int> GetIndexMatrix() const
      {
        return m_IndexMatrixOnDevice;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief sorts k values of columns in matrix into first k positions
      //!        also provides buffer of keys as "index buffer" through m_IndexBuffer variable which can be retrieved
      //!        using Get function
      //! @param DATA matrix buffer to sort
      //! @param NR_TO_SORT the number of values to sort
      //! @return buffer with first k lowest values sorted in first k columns
      Matrix< float> operator()
      (
        Matrix< float> &DATA, const size_t &NR_TO_SORT
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class InsertionSort

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_INSERTION_SORT_H_ 
