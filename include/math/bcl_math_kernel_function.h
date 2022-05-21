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

#ifndef BCL_MATH_KERNEL_FUNCTION_H_
#define BCL_MATH_KERNEL_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KernelFunction
    //! @brief represents a function that is 1 between 0 and 1, 0 otherwise
    //!
    //! @see @link example_math_kernel_function.cpp @endlink
    //! @author mueller
    //! @date Apr 26, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KernelFunction :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      //! range of the kernel
      Range< double> m_Kernel;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      KernelFunction();

      //! @brief construct from properties
      KernelFunction( const Range< double> &KERNEL);

      //! @brief Clone function
      //! @return pointer to new KernelFunction
      KernelFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief access to range of the kernel
      //! @return the range of the kernel
      Range< double> GetKernel() const;

      //! @brief access to range of the kernel
      //! @param KERNEL the new range of the kernel
      void SetKernel( const Range< double> &KERNEL);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator f( x) = y = 1 if x is in the kernel range
      //! @param ARGUMENT the x
      //! @return f( x)
      double operator()( const double &ARGUMENT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class KernelFunction

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_KERNEL_FUNCTION_H_ 
