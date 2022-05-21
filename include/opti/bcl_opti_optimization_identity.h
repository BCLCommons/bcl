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

#ifndef BCL_OPTI_OPTIMIZATION_IDENTITY_H_
#define BCL_OPTI_OPTIMIZATION_IDENTITY_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_optimization_interface.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OptimizationIdentity
    //! @brief Class for testing that does not alter the argument
    //!
    //! @see @link example_opti_optimization_identity.cpp @endlink
    //! @author fischea
    //! @date Nov 08, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ArgumentType>
    class OptimizationIdentity :
      public OptimizationInterface< t_ArgumentType>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      OptimizationIdentity()
      {
      }

      //! @brief clone function
      //! @return pointer to a new OptimizationIdentity
      OptimizationIdentity *Clone() const
      {
        return new OptimizationIdentity( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "OptimizationIdentity");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Implementation for testing, which does not alter the provided argument.");

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief optimizes the provided argument
      //! @detail the argument is not altered
      //! @param ARGUMENT data to be optimized
      void Optimize( t_ArgumentType &ARGUMENT) const
      {
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class OptimizationIdentity

    //! single instance of this class
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> OptimizationIdentity< t_ArgumentType>::s_Instance
    (
      util::Enumerated< OptimizationInterface< t_ArgumentType> >::AddInstance( new OptimizationIdentity< t_ArgumentType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_OPTIMIZATION_IDENTITY_H_
