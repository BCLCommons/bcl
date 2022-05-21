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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math_kernel_function.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> KernelFunction::s_Instance
    (
      GetObjectInstances().AddInstance( new KernelFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KernelFunction::KernelFunction()
    {
    }

    //! @brief construct from properties
    KernelFunction::KernelFunction
    (
      const Range< double> &KERNEL
    ) :
      m_Kernel( KERNEL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new KernelFunction
    KernelFunction *KernelFunction::Clone() const
    {
      return new KernelFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &KernelFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to range of kernel
    //! @return the range of the kernel
    Range< double> KernelFunction::GetKernel() const
    {
      return m_Kernel;
    }

    //! @brief access to range of kernel
    //! @param KERNEL the new range of the kernel
    void KernelFunction::SetKernel( const Range< double> &KERNEL)
    {
      m_Kernel = KERNEL;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator f( x) = y = 1 if x is in the kernel range
    //! @param ARGUMENT the x
    //! @return f( x)
    double KernelFunction::operator()( const double &ARGUMENT) const
    {
      return ( m_Kernel.IsWithin( ARGUMENT)) ? 1.0 : 0.0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelFunction::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &KernelFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
