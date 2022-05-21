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

#ifndef BCL_CHEMISTRY_ROTAMER_CLUSTER_CENTER_H_
#define BCL_CHEMISTRY_ROTAMER_CLUSTER_CENTER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_priority_dihedral_angles.h"
#include "bcl_chemistry_rotamer_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotamerClusterCenter
    //! @brief Class that finds conformer/rotamer cluster centers for all unique conformations of a given configuration
    //!
    //! @see @link example_chemistry_rotamer_cluster_center.cpp @endlink
    //! @author kothiwsk
    //! @date Aug 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotamerClusterCenter :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      //! configuration of interest whose rotamer library has to be built
      util::ShPtr< FragmentConfigurationShared>     m_Molecule;

      //! conformation comparison object
      ConformationComparisonByDihedralBins          m_ConformationComparison;

      //! priority dihedral angle object of configuration of interest
      PriorityDihedralAngles                        m_PriorityDihedral;

      //! storage container to store different conformation clusters
      storage::List< RotamerEnsemble>               m_Clusters;

      //! bin size used for conformation comparison
      double                                        m_BinSize;

      //! Number of states that were never seen but were simulated
      double                                        m_NumberSimulatedButNotSeen;

      //! Whether this is a pure ring system
      bool                                          m_IsPureRing;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param MOLECULE the configuration of interest
      //! @param TOLERANCE the tolerance level to be used for dihedral binning
      RotamerClusterCenter( const ConformationInterface &MOLECULE, const double BIN_SIZE = double( 30.0));

      //! @brief Clone function
      //! @return pointer to new RotamerClusterCenter
      RotamerClusterCenter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the configuration whose conformation clusters this object calculates
      //! @return configuration whose conformation clusters this object calculates
      const FragmentConfigurationShared &GetConfiguration() const;

      //! @brief returns reference to all conformation clusters
      //! @return reference to all conformation clusters
      math::RunningAverageSD< double> GetStats() const;

      //! @brief get the count of rotamer seen most number of times
      //! @return count of rotamer seen most number of times
      double GetMaxWeight() const;

      //! @brief get the count of rotamer seen most number of times
      //! @return count of rotamer seen most number of times
      double GetSimulatedWeightUnseenRotamers() const;

      //! @brief get the count of rotamer seen most number of times
      //! @return count of rotamer seen most number of times
      size_t GetMaxInstance() const;

      //! @brief returns reference to all conformation clusters
      //! @return reference to all conformation clusters
      const storage::List< RotamerEnsemble> &GetClusters() const;

      //! @brief returns reference to all conformation clusters
      //! @return reference to all conformation clusters
      const PriorityDihedralAngles &GetPriority() const
      {
        return m_PriorityDihedral;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculates which bin signature the given conformation belongs to
      //! @param CONFORMATION conformation to be used for calculation
      //! @return void
      void operator()
      (
        const FragmentComplete &CONFORMATION,
        const double WEIGHT,
        const bool &SIMULATED = false
      );

      //! @brief calculates which bin signature the given conformation belongs to
      //! @param CONFORMATION conformation to be used for calculation
      //! @return void
      void operator()
      (
        const graph::SubgraphIsomorphism< util::SiPtr< const AtomConformationalInterface>, size_t> &ATOM_ISO,
        const storage::Vector< storage::Vector< size_t> > &ISO,
        const std::string &NAME,
        const bool &SIMULATED = false
      );

      //! @brief calculates the centers of cluster from ensembles stored in m_ConformationEnsembles
      //! @return void
      void CalculateClusterCenters( size_t MAX_COUNTS);

      //! @brief get an average structure to get a good structure
      //! @return void
      FragmentConformationShared GetAverageStucture( size_t MAX_COUNTS);

      //! @brief removes those rotamers which whose weight is less than cut off value
      //! @return void
      void PruneRotamersByWeight( double CUT_OFF);

      void PruneRotamersByInstances( size_t CUT_OFF);

      //! @brief removes those rotamers which have non-planar aromatic rings
      //! @return void
      void PruneNonPlanarAromaticRings();

      //! @brief get weights associated with each isomorphism
      static storage::Vector< double> GetIsomorphismWeights
      (
        const storage::Vector< storage::Vector< size_t> > &ISO,
        const size_t &MOLECULE_SIZE
      );

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

    private:

      //! @brief Create a graph of a conformation
      //! @param CONFORMATION a conformation of configuration of interest
      //! @param DIHEDRAL_BIN dihedral bin signature of the conformation being passed in
      //! @param ATOM_DIHEDRAL_INDICES the atoms that make the priority dihedral angles for the conformation being passed in
      //! @return The conformation converted into a graph with coloring scheme that takes into account dihedral bin
      graph::ConstGraph< size_t, size_t> CreateGraphWithDihedralBinsOnEdges
      (
        const ConformationInterface &CONFORMATION,
        const storage::Vector< int> &DIHEDRAL_BINS,
        const storage::Vector< storage::VectorND< 4, size_t> > &ATOM_DIHEDRAL_INDICES
      ) const;

    }; // class RotamerClusterCenter
  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ROTAMER_CLUSTER_CENTER_H_
