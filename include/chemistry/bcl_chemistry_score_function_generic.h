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

#ifndef BCL_CHEMISTRY_SCORE_FUNCTION_GENERIC_H_
#define BCL_CHEMISTRY_SCORE_FUNCTION_GENERIC_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_conformation_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreFunctionGeneric
    //! @brief Initializes a score function as a CheminfoProperty and returns the mean property value across all indices
    //!
    //! @see @link example_chemistry_score_function_generic.cpp @endlink
    //! @author brownbp1
    //! @date May 01, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreFunctionGeneric :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {
    private:

      // the descriptor to use
      descriptor::CheminfoProperty m_Descriptor;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ScoreFunctionGeneric();

      //! @brief constructor with parameters
      //! @param DESCRIPTOR the descriptor to use
      ScoreFunctionGeneric
      (
        const descriptor::CheminfoProperty &DESCRIPTOR
      );

      //! @brief Clone function
      //! @return pointer to new ScoreFunctionGeneric
      ScoreFunctionGeneric *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the class name when used in a dynamic context
      //! @return the class name when used in a dynamic context
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return propensity score for observing rotamers that exist in conformation
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief gets a serializer for constructing this class in a dynamic context
      //! @return the serializer containing member data
      io::Serializer GetSerializer() const;

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

    }; // class ScoreFunctionGeneric

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SCORE_FUNCTION_GENERIC_H_
