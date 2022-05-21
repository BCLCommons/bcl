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

#ifndef BCL_CHEMISTRY_ROTAMER_ENSEMBLE_H_
#define BCL_CHEMISTRY_ROTAMER_ENSEMBLE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_conformation_shared.h"
#include "bcl_chemistry_priority_dihedral_angles.h"
#include "graph/bcl_graph_const_graph.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotamerEnsemble
    //! @brief Class that stores all conformations that belong to one cluster or have the same dihedral signature
    //!
    //! @see @link example_chemistry_rotamer_ensemble.cpp @endlink
    //! @author kothiwsk
    //! @date Aug 07, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotamerEnsemble :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      //! Storage container to store all conformations having same dihedral bin signature
      storage::List< FragmentConformationShared>        m_Conformers;

      //! running average and standard deviation to collect statistics about cluster of conformations
      math::RunningAverage< linal::Vector< double> > m_AverageBinAngles;

      //! graph colored by dihedral bins
      graph::ConstGraph< size_t, size_t>                m_Graph;

      //! cluster center
      FragmentConformationShared                        m_Center;

      //! bin size used for conformation comparison
      double                                            m_BinSize;

      //! get number of instances of the conformation
      size_t                                            m_Instances;

      //! # of simulated instances
      double                                            m_SimulatedInstances;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      RotamerEnsemble()
      {
      }

      //! @brief constructor
      //! @param MOLECULE Molecule that needs to be added to the cluster
      //! @param GRAPH method to be used for finding cluster center
      //! @param BIN_SIZE dihedral bin size used for conformation determination
      RotamerEnsemble
      (
        const FragmentConformationShared &MOLECULE,
        const graph::ConstGraph< size_t, size_t> &GRAPH,
        const double WEIGHT,
        double BIN_SIZE,
        const bool &SIMULATED
      );

      //! @brief Clone function
      //! @return pointer to new RotamerEnsemble
      RotamerEnsemble *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns average dihedral angles for conformations that belong to this cluster
      //! @return average dihedral angles
      const linal::Vector< double> &GetAverage() const;

      //! @brief returns the size of this cluster
      //! @return the size of the cluster
      const double &GetWeights() const;

      //! @brief returns the number of times fragment has been seen in the structure database
      //! @return the number of times fragment has been seen in the structure database
      const size_t &GetNumberInstances() const;

      //! @brief return the conformation graph colored by dihedral bins for this cluster
      //! @return the conformatin graph for this cluster
      const graph::ConstGraph< size_t, size_t> &GetGraph() const;

      //! @brief return the conformation which is the center of this cluster
      //! @return the center of this cluster
      const FragmentConformationShared &GetClusterCenter() const;

      //! @brief return the conformation which is the center of this cluster
      //! @return the center of this cluster
      const double &GetSimulatedInstances() const
      {
        return m_SimulatedInstances;
      }

      //! @brief return members of the rotamer ensemble
      //! @return all members of rotamer ensemble
      const storage::List< FragmentConformationShared> &GetRotamers() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the center of this cluster from the cluster info
      //! @param MAX_COUNTS maximum number of examples over which to calculate cluster center
      //! @return void
      void CalculateCenter( size_t MAX_COUNTS);

    ///////////////
    // operators //
    ///////////////

      //! @brief add a conformation to the cluster info
      //! @param FRAGMENT the conformation to add
      void operator()( const FragmentConformationShared &FRAGMENT, const double WEIGHT, const bool &SIMULATED = false);

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
    };
  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ROTAMER_ENSEMBLE_H_
