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

#ifndef BCL_MATH_MUTATE_INTERFACE_H_
#define BCL_MATH_MUTATE_INTERFACE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateInterface
    //! @brief interface for a Monte Carlo based Minimization on t_ResultType = Function( t_ArgumentType)
    //! @details it provides an t_ArgumentType operator()( t_ArgumentType) that alters the ARGUMENT.
    //!
    //! @remarks example unnecessary
    //! @author woetzen, karakam
    //! @date 01.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class MutateInterface :
      public FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> >
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      virtual ~MutateInterface()
      {
      }

      //! @brief clone function
      //! @return pointer to a new MutateInterface
      virtual MutateInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const = 0;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      virtual const std::string &GetScheme() const
      {
        return GetClassIdentifier();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param ARGUMENT Argument of interest
      //! @return MutateResult that results from mutating to the argument
      virtual MutateResult< t_ArgumentType> operator()( const t_ArgumentType &ARGUMENT) const = 0;

    }; // template class MutateInterface

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_MUTATE_INTERFACE_H_
