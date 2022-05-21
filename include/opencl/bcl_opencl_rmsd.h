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

#ifndef BCL_OPENCL_RMSD_H_
#define BCL_OPENCL_RMSD_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_vector.h"
#include "linal/bcl_linal_matrix.h"
#include "model/bcl_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RMSD
    //! @brief performs rmsd calculation given two arrays
    //!
    //! @see @link example_opencl_rmsd.cpp @endlink
    //! @author loweew
    //! @date Mar 23, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RMSD :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

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

      //! @brief default constructor
      RMSD();

      //! @brief constructor from command queue
      //! @param QUEUE command queue
      RMSD( const CommandQueue &QUEUE);

      //! @brief Clone function
      //! @return pointer to new RMSD
      RMSD *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief rmsd calculation
      //! @param INPUT_A, INPUT_B cl buffers of matrices to compare
      //! @return rmsd of the INPUT_A vs INPUT_B
      float operator()
      (
        const Matrix< float> &INPUT_A,
        const Matrix< float> &INPUT_B
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

    }; // class RMSD

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_RMSD_H_ 
