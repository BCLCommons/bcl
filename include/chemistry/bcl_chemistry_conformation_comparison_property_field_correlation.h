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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_PROPERTY_FIELD_CORRELATION_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_PROPERTY_FIELD_CORRELATION_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "util/bcl_util_implementation.h"
namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonPropertyFieldCorrelation
    //! @brief This class is designed to be used for determining and comparing 3D structures for molecules based
    //!        on euclidean and property distance, and even for molecules with differing constitutions
    //!
    //! @see @link example_chemistry_conformation_comparison_property_field_correlation.cpp @endlink
    //! @author mendenjl, brownbp1
    //! @date Nov 15, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonPropertyFieldCorrelation :
      public ConformationComparisonInterface
    {

    protected:
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! Atom properties to include in the distance
      mutable descriptor::Combine< AtomConformationalInterface, float> m_Properties;

      //! Weights for each atom property in the overall correlation
      linal::Vector< float> m_PropertyWeights;
      linal::Vector< float> m_PropertyWeightsSqr;

      //! max distance two atoms on different molecules can be apart and still considered for property distance computations
      double m_MaxAtomDistanceCriterion;

      //! native poses
      mutable util::SiPtr< const FragmentEnsemble> m_EnsembleB;

      //! whether to perform molecule:Compare w/ SymmetryRMSD against native poses in EnsembleB
      bool m_OptimizingWeights;

      //! for when m_OptimizingWeights is set
      mutable storage::Vector< ConformationComparisonBySymmetryRmsd> m_SymmetryRmsdComparers;

      //! path for .csv output for when m_OptimizingWeights is set
      std::string m_Path;

      //! filename to be appended to the path for when m_OptimizingWeights is set
      std::string m_FileName;

      //! Linear penalty for unmatched atoms
      double m_MismatchPenalty;

      //! Fraction below which the heavy mismatch penalty is applied
      double m_HeavyPenaltyFraction;

      //! Penalty for having < m_HeavyPenaltyFraction of atoms mutually matched on the maximally matched molecule
      double m_HeavyMismatchPenalty;

      // TODO: I do not believe we need/use this anymore for this score
      //! Anchor weight to avoid molecules floating away
      double m_AnchorWeight;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConformationComparisonPropertyFieldCorrelation();

      //! @brief full constructor
      ConformationComparisonPropertyFieldCorrelation
      (
        const descriptor::Combine< AtomConformationalInterface, float> &PROPERTIES,
        const linal::Vector< float> &PROPERTY_WEIGHTS,
        const double MAX_ATOM_DIST,
        const bool OPTI_WEIGHTS,
        const double MISMATCH_PENALTY,
        const double HEAVY_MISMATCH_FRACTION,
        const double HEAVY_MISMATCH_PENALTY,
        const double ANCHOR_WEIGHT
      );

      //! virtual copy constructor
      ConformationComparisonPropertyFieldCorrelation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the max ratio of matched atoms between molecules a and b
      const double GetMaxAtomDistance() const;

      //! @brief get the linear penalty for unmatched atoms
      const double GetLinearPenalty() const;

      //! @brief get the fraction below which the heavy mismatch penalty is applied
      const double GetHeavyPenaltyFraction() const;

      //! @brief get penalty for having < m_HeavyPenaltyFraction of atoms mutually matched on the maximally matched molecule
      const double GetHeavyPenalty() const;

//      //! @brief get the anchor weight to avoid molecules floating away
//      const double GetAnchorWeight() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief find the RMSDX between two conformations
      //! @param MOLECULE_A - the fragment to align against
      //! @param MOLECULE_B - the fragment being aligned
      //! @return the RMSDX between MOLECULE_A and MOLECULE_B
      double operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

      //! @brief find the RMSDX between two conformations
      //! @param MOLECULE_A - the fragment to align against
      //! @param MOLECULE_B - the fragment being aligned
      //! @param MAX_ATOM_DISTANCE - the maximum atom distance determining atom pairs
      //! @return the RMSDX between MOLECULE_A and MOLECULE_B
      storage::Pair< double, double> PropertyDistanceScorer
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B,
        const float &MAX_ATOM_DISTANCE,
        const storage::Vector< size_t> &EXCLUDE_ATOMS_A = storage::Vector< size_t>(),
        const storage::Vector< size_t> &EXCLUDE_ATOMS_B = storage::Vector< size_t>()
      ) const;

      //! @brief obtain mutually closest atom pairs of two molecules
      //! @param MOLECULE_A - first input molecule
      //! @param MOLECULE_B - second input molecule
      //! @param MOL_A_ALIGNED_ATOMS - indices of aligned atoms in molecule A
      //! @param MOL_B_ALIGNED_ATOMS - indices of aligned atoms in molecule B
      //! @return void
      static storage::Pair< storage::Vector< size_t>, storage::Vector< size_t> > GetAlignedAtoms
      (
        const FragmentComplete &MOLECULE_A,
        const FragmentComplete &MOLECULE_B,
        const storage::Vector< size_t> &KEEP_INDICES_A,
        const storage::Vector< size_t> &KEEP_INDICES_B
      );

      //! @brief find parameters for given atom coordinates
      //! @param ATOM_COORD_X - the x coordinate
      //! @param ATOM_COORD_Y - the y coordinate
      //! @param ATOM_COORD_Z - the z coordinate
      //! @return the values of properties at each 3D coordinate
      linal::Vector< float> GetPropertyAtCoordinates
      (
        const ConformationInterface &MOLECULE
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief prepare the class for comparing a conformation
      //! @param MOLECULE the molecule to prepare to compare
      void Prepare( const ConformationInterface &MOLECULE) const;

      //! @brief prepare for comparison against an ensemble
      //! @param ENSEMBLE the ensemble against which MOLECULE will be compared
      void PrepareEnsemble( const FragmentEnsemble &ENSEMBLE) const;

      //! @brief set atom properties and remove size zero properties
      //! @param PROPERTIES descriptor label to use
      descriptor::Combine< AtomConformationalInterface, float> PrepareProperties( const util::ObjectDataLabel &PROPERTIES) const;

      //! @brief sets properties with which to compare alignment poses
      void SetProperties( const descriptor::Combine< AtomConformationalInterface, float> &PROPERTIES);

      //! @brief sets the weights of each property
      void SetPropertyWeights( const linal::Vector< float> &PROPERTY_WEIGHTS);

      //! @brief get the max ratio of matched atoms between molecules a and b
      void SetMaxAtomDistance( const double MAX_ATOM_DISTANCE);

      //! @brief get the linear penalty for unmatched atoms
      void SetLinearPenalty( const double MISMATCH_PENALTY);

      //! @brief get the fraction below which the heavy mismatch penalty is applied
      void SetHeavyPenaltyFraction( const double HEAVY_PENALTY_FRACTION);

      //! @brief get penalty for having < m_HeavyPenaltyFraction of atoms mutually matched on the maximally matched molecule
      void SetHeavyPenalty( const double HEAVY_MISMATCH_PENALTY);

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );
    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_PROPERTY_FIELD_CORRELATION_H_
