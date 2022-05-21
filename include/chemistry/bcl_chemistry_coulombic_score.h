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

#ifndef BCL_CHEMISTRY_COULOMBIC_SCORE_H_
#define BCL_CHEMISTRY_COULOMBIC_SCORE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "descriptor/bcl_descriptor_base_sequence.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CoulombicScore
    //! @brief Scores the conformation based on propensity of observing rotamers of constituent fragments
    //! @details for a given conformation, score is calculated based on which rotamers have been used for sampling
    //!           conformations
    //!
    //! @see @link example_chemistry_coulombic_score.cpp @endlink
    //! @author kothiwsk
    //! @date Jan 05, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CoulombicScore :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {
    private:

        descriptor::CheminfoProperty                                          m_Charge;

        double                                                                m_Weight;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      CoulombicScore()
      {
      }

      CoulombicScore( const descriptor::CheminfoProperty &ATOM_CHARGE, double WEIGHT);

      //! @brief Clone function
      //! @return pointer to new CoulombicScore
      CoulombicScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

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

    //////////////////////
    // helper functions //
    //////////////////////
    }; // class CoulombicScore

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_COULOMBIC_SCORE_H_
