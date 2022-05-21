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

#ifndef BCL_CHEMISTRY_ATOM_ONE_FOUR_INTERACTION_SCORE_H_
#define BCL_CHEMISTRY_ATOM_ONE_FOUR_INTERACTION_SCORE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomOneFourInteractionScore
    //! @brief Scores clashes between atoms of a given conformation
    //! @details Potential that evalutes clashes in a given molecule conformation
    //!
    //! @see @link example_chemistry_atom_one_four_interaction_score.cpp @endlink
    //! @author mendenjl
    //! @date Sep 26, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomOneFourInteractionScore :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {

    private:

      //! Weight for length for this score.
      double m_LengthWeight;

      //! Weight for identity of interactions
      double m_InteractionIdentityWeight;

      //! Molecule dependent variables. These variables are cached because they are expensive to compute and the same
      //! molecule may have the clash score called on it hundreds of times in succession
      mutable storage::Vector< AtomType> m_LastAtomTypes;
      mutable storage::Vector< sdf::BondInfo> m_LastBondInfo;

      //! Current scores
      mutable storage::Vector< double> m_AtomOneFourInteractionScores;

      // 1-4 neighbors that are not part of the same rigid component
      mutable storage::Vector
      <
        storage::Vector< storage::Pair< size_t, util::SiPtr< const storage::VectorND< 3, double> > > >
      > m_OneFourNeighbors;

      //! @brief get the map from Key to average interaction frequency, average distance, and standard deviation of distance
      static storage::Map< std::string, storage::VectorND< 3, double> > s_InteractionMap;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief default constructor. Can take whether to consider hydrogens or tolerance
      //! @param HTOL hydrogen VDW-overlap tolerance. Default value causes non-covalent overlap of hydrogens to be ignored
      explicit AtomOneFourInteractionScore( double LENGTH_WEIGHT = 1.0, double INTERACTION_WEIGHT = 1.0);

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AtomOneFourInteractionScore
      AtomOneFourInteractionScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return clash score for the given atom pair
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief evaluate clashes for given atoms
      //! @param MOLECULE molecule that needs to scored
      //! @return clash score for the given atoms
      double operator()
      (
        const ConformationInterface &MOLECULE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief estimate the minimum distance of two atoms of a given type, bond distance, and presence in the same rigid structure
      static double EstimateMinimumDistance
      (
        const AtomConformationalInterface &A,
        const AtomConformationalInterface &B,
        const int &BOND_DISTANCE,
        const storage::Vector< size_t> &RIGID_SUBSTRUCTURES_A,
        const storage::Vector< size_t> &RIGID_SUBSTRUCTURES_B
      );

      //! @brief initialize the interaction map
      //! this is slow and should be performed only once
      static bool InitializeInteractionMap();

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief Test whether this molecule is the same (constitutionally) as the molecule for which the state of this
      //!        class currently can handle
      bool IsMoleculeInfoCached( const ConformationInterface &CONF) const;

      //! @brief Update molecule change the molecule that this class will compute the clash score for
      //! @param MOL molecule of interest
      void UpdateMolecule( const ConformationInterface &CONF) const;

    }; // class AtomOneFourInteractionScore

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_ONE_FOUR_INTERACTION_SCORE_H_
