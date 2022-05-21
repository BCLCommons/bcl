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
#include "chemistry/bcl_chemistry_conformation_comparison_by_symmetry_rmsd.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

//#define BCL_PROFILE_SymmetryRMSD
#ifdef BCL_PROFILE_SymmetryRMSD
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonBySymmetryRmsd::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonBySymmetryRmsd())
    );
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonBySymmetryRmsd::s_RealSpaceInstance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonBySymmetryRmsd( false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, accepts bool that indicates whether to superimpose (default) or not
    ConformationComparisonBySymmetryRmsd::ConformationComparisonBySymmetryRmsd( const bool &SUPERIMPOSE, const size_t &ISO_LIMIT) :
      m_Superimpose( SUPERIMPOSE),
      m_ConsiderH( false),
      m_AtomComparison( ConformationGraphConverter::e_ElementType),
      m_BondComparison( ConfigurationalBondTypeData::e_FuzzyBondOrderAmideOrAromaticWithRingness),
      m_NumberIsomorphismLimit( ISO_LIMIT)
    {
    }

    //! virtual copy constructor
    ConformationComparisonBySymmetryRmsd *ConformationComparisonBySymmetryRmsd::Clone() const
    {
      return new ConformationComparisonBySymmetryRmsd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonBySymmetryRmsd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonBySymmetryRmsd::GetAlias() const
    {
      static const std::string s_Name( "SymmetryRMSD"), s_real_name( "SymmetryRealSpaceRMSD");
      return m_Superimpose ? s_Name : s_real_name;
    }

    //! @brief returns the number of isomorphisms
    //! @return the number of isomorphisms
    size_t ConformationComparisonBySymmetryRmsd::GetNumberIsomorphisms() const
    {
      return m_Isomorphisms.GetSize();
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonBySymmetryRmsd::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      // Determine if conformations can be compared
      if( !ConformationComparisonInterface::ConformationsAreComparable( MOLECULE_A, MOLECULE_B))
      {
        // nope, so return an undefined
        return util::GetUndefined< double>();
      }

      // if there were 0-1 atoms in the molecules, then the rmsd is automatically 0
      if( m_Superimpose && MOLECULE_A.GetNumberAtoms() < size_t( 2))
      {
        return 0.0;
      }

      if
      (
        m_Isomorphisms.IsEmpty()
        || m_PreviousAtomTypes != MOLECULE_A.GetAtomTypesVector()
        || m_PreviousEdgeTypes != MOLECULE_A.GetAdjacencyList( m_BondComparison)
      )
      {
        Prepare( MOLECULE_A);
      }

      if
      (
        MOLECULE_A.GetChiralityString() + " " + MOLECULE_A.GetDoubleBondIsometryString()
        != MOLECULE_B.GetChiralityString() + " " + MOLECULE_B.GetDoubleBondIsometryString()
      )
      {
        BCL_MessageStd
        (
          "Comparing molecules with differing chirality or isometry: " +
          MOLECULE_A.GetChiralityString() + " " + MOLECULE_A.GetDoubleBondIsometryString()
          + " vs " + MOLECULE_B.GetChiralityString() + " " + MOLECULE_B.GetDoubleBondIsometryString()
        );
      }

      const storage::Vector< quality::RMSDPreprocessor> &a_data( GetCachedInfo( MOLECULE_A));

      // if there were 0-1 atoms in the molecules, then the rmsd is automatically 0
      if( a_data( 0).GetSize() < size_t( 2))
      {
        return 0.0;
      }

      const quality::RMSDPreprocessor &b_preproc( GetCachedInfo( MOLECULE_B)( 0));

      const double min_rmsd( b_preproc.GetSize() * std::numeric_limits< double>::epsilon());
      double rmsd( math::GetHighestBoundedValue< double>());
      for
      (
        auto itr_iso( a_data.Begin()), itr_iso_end( a_data.End());
        itr_iso != itr_iso_end;
        ++itr_iso
      )
      {

        // find the RMSDs between the coordinate vectors
        double cur_rmsd
        (
          m_Superimpose
          ? itr_iso->SuperimposedRMSD( b_preproc)
          : itr_iso->RMSD( b_preproc)
        );

        if( util::IsDefined( cur_rmsd) && cur_rmsd < rmsd)
        {
          rmsd = cur_rmsd;
          if( rmsd <= min_rmsd)
          {
            break;
          }
        }
      }

      // if the rmsd is very close to zero, set it to zero
      // assume that errors of size ~std::numeric_limits< double>::epsilon() are accumulated in the RMSD calculation for every atom
      if( math::EqualWithinAbsoluteTolerance( 0.0, rmsd, b_preproc.GetSize() * std::numeric_limits< double>::epsilon()))
      {
        rmsd = 0.0;
      }

      return rmsd;
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void ConformationComparisonBySymmetryRmsd::Prepare( const ConformationInterface &MOLECULE) const
    {
      #ifdef BCL_PROFILE_SymmetryRMSD
      static util::Stopwatch s_timer( "symmetry rmsd::prepare", util::Message::e_Standard, true);
      s_timer.Start();
      #endif
      graph::ConstGraph< size_t, size_t> molecule_graph
      (
        ConformationGraphConverter
        (
          m_AtomComparison,
          m_BondComparison,
          !m_ConsiderH
        )( MOLECULE)
      );
      graph::SubgraphIsomorphism< size_t, size_t> csi_substructure;
      csi_substructure.SetGraphs( molecule_graph, molecule_graph);

      csi_substructure.FindAllIsomorphisms();

      m_Isomorphisms = csi_substructure.GetIsomorphisms();
      BCL_MessageStd( "Found " + util::Format()( m_Isomorphisms.GetSize()) + " isomorphisms for this molecule");
      if( m_Isomorphisms.GetSize() > m_NumberIsomorphismLimit)
      {
        std::random_shuffle( m_Isomorphisms.Begin() + 1, m_Isomorphisms.End(), random::GetRandomSizeT);
        m_Isomorphisms.Resize( m_NumberIsomorphismLimit);
      }
      m_PreviousAtomTypes = MOLECULE.GetAtomTypesVector();
      m_PreviousEdgeTypes = MOLECULE.GetAdjacencyList( m_BondComparison);
      #ifdef BCL_PROFILE_SymmetryRMSD
      s_timer.Stop();
    #endif
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonBySymmetryRmsd::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the root-mean-squared-deviation between two molecules that are identical on"
        "the constitutional level, and whose atoms have already been aligned"
      );
      parameters.AddInitializer
      (
        "consider H",
        "set to true if you want to consider RMSD of hydrogens atoms as well. This is usually a bad idea in larger molecules, "
        "because the number of possible automorphisms of the molecule with hydrogens grows exponentially. "
        "Only recommended for molecules with fewer than 15 heavy atoms",
        io::Serialization::GetAgent( &m_ConsiderH),
        "false"
      );
      parameters.AddInitializer
      (
        "atom comparison",
        "Used to decide which atom types can be considered equivalent",
        io::Serialization::GetAgent( &m_AtomComparison),
        ConformationGraphConverter::GetAtomComparisonType( ConformationGraphConverter::e_ElementType)
      );
      parameters.AddInitializer
      (
        "bond comparison",
        "Used to decide whether bonds are similar enough to consider equivalent",
        io::Serialization::GetAgent( &m_BondComparison),
        ConfigurationalBondTypeData::GetDataName( ConfigurationalBondTypeData::e_FuzzyBondOrderAmideOrAromaticWithRingness)
      );
      parameters.AddInitializer
      (
        "isomorphism limit",
        "If the calculation is prohibitively slow for highly-symmetric molecules, set the isomorphism limit to speed up the "
        "computation. There are exponentially diminishing returns for having all isomorphisms present, so usually 10-100 is "
        "sufficient for most purposes",
        io::Serialization::GetAgent( &m_NumberIsomorphismLimit),
        util::Format()( 1000)
      );
      return parameters;
    }

    //! @brief prepare the cache entry for a given molecule
    const storage::Vector< quality::RMSDPreprocessor> &ConformationComparisonBySymmetryRmsd::GetCachedInfo( const ConformationInterface &MOLECULE) const
    {
      auto itr( m_Cache.Find( size_t( &MOLECULE)));
      if( itr != m_Cache.End())
      {
        return itr->second;
      }
      storage::Vector< quality::RMSDPreprocessor> &new_vec
      (
        m_Cache.Insert
        (
          std::make_pair
          (
            size_t( &MOLECULE),
            storage::Vector< quality::RMSDPreprocessor>()
          )
        ).first->second
      );
      new_vec.AllocateMemory( m_Isomorphisms.GetSize());
      auto initial_coords( m_ConsiderH ? MOLECULE.GetAtomCoordinates() : MOLECULE.GetHeavyAtomCoordinates());
      for( auto itr_iso( m_Isomorphisms.Begin()), itr_iso_end( m_Isomorphisms.End()); itr_iso != itr_iso_end; ++itr_iso)
      {
        auto initial_coords_copy( initial_coords);
        initial_coords_copy.Reorder( *itr_iso);
        new_vec.PushBack( quality::RMSDPreprocessor( initial_coords_copy, m_Superimpose));
      }
      MOLECULE.GetChangeSignal().Connect
      (
        const_cast< ConformationComparisonBySymmetryRmsd *>( this),
        &ConformationComparisonBySymmetryRmsd::RemoveResultFromCache
      );
      return new_vec;
    }

    //! @brief remove results for this object from the cache
    //! @param ARGUMENT argument that should be removed
    void ConformationComparisonBySymmetryRmsd::RemoveResultFromCache( const ConformationInterface &ARGUMENT)
    {
      // find result for that Argument
      auto itr( m_Cache.Find( size_t( &ARGUMENT)));

      // skip if there is no result cached
      if( itr != m_Cache.End())
      {
        m_Cache.RemoveElement( itr);
        ARGUMENT.GetChangeSignal().Disconnect( const_cast< ConformationComparisonBySymmetryRmsd *>( this));
      }
      else
      {
        BCL_MessageStd( "hunh; this shouldn't happen");
      }
    }

  } // namespace chemistry
} // namespace bcl
