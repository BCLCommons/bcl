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
#include "chemistry/bcl_chemistry_minimize_ligand_in_pocket.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_clash_score.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_ligand_pocket_fit_score.h"
#include "chemistry/bcl_chemistry_mutate_molecule_generic.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MinimizeLigandInPocket::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MinimizeLigandInPocket())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MinimizeLigandInPocket::MinimizeLigandInPocket() :
      ConformationComparisonPropertyFieldCorrelation(),
      m_Iterations( 500),
      m_IterationLimit( 40),
      m_StartPosition( util::GetUndefinedDouble()),
      m_Movie( "")
    {
    }

    //! @brief full constructor
    MinimizeLigandInPocket::MinimizeLigandInPocket
    (
      const size_t &ITERATIONS,
      const size_t &ITERATION_LIMIT,
      const linal::Vector3D &START_POSITION
    ) :
      ConformationComparisonPropertyFieldCorrelation(),
      m_Iterations( ITERATIONS),
      m_IterationLimit( ITERATION_LIMIT),
      m_StartPosition( START_POSITION),
      m_Movie( "")
    {
    }

    //! @brief movie constructor
    MinimizeLigandInPocket::MinimizeLigandInPocket
    (
      const size_t &ITERATIONS,
      const size_t &ITERATION_LIMIT,
      const linal::Vector3D &START_POSITION,
      const std::string &MOVIE_FILENAME
    ) :
      ConformationComparisonPropertyFieldCorrelation(),
      m_Iterations( ITERATIONS),
      m_IterationLimit( ITERATION_LIMIT),
      m_StartPosition( START_POSITION),
      m_Movie( MOVIE_FILENAME)
    {
    }

    //! virtual copy constructor
    MinimizeLigandInPocket *MinimizeLigandInPocket::Clone() const
    {
      return new MinimizeLigandInPocket( *this);
    }

  /////////////////
  // data access // 
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MinimizeLigandInPocket::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &MinimizeLigandInPocket::GetAlias() const
    {
      static std::string s_name( "LigFit");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the property RMSD
    //! @param MOLECULE - molecule to be fit into pocket
    //! @param POCKET - pocket into which MOLECULE is being geometrically fit
    //! @return the molecule in its new pose relative to the static pocket
    double MinimizeLigandInPocket::operator()
    (
      const ConformationInterface &MOLECULE,
      const ConformationInterface &POCKET
    ) const
    {
      storage::Pair< FragmentComplete, double> ligand_pocket_orientation_score;
      if( m_StartPosition.IsDefined())
      {
        storage::Pair< FragmentComplete, double> ligand_pocket_orientation_score( MinimizeLigandInPocket::OrientLigandInPocketEngine( MOLECULE, POCKET, m_Iterations, m_IterationLimit, m_StartPosition));
      }
      else
      {
        storage::Pair< FragmentComplete, double> ligand_pocket_orientation_score( MinimizeLigandInPocket::OrientLigandInPocketEngine( MOLECULE, POCKET, m_Iterations, m_IterationLimit, POCKET.GetCenter()));
      }

      //output best molecule poses with scores
      io::OFStream outfile;
      io::File::MustOpenOFStream( outfile, "minimize_ligand_in_pocket_pilot.sdf", std::ios::app);
      ligand_pocket_orientation_score.First().StoreProperty( "LigPocketCollisionScore", util::Format()( ligand_pocket_orientation_score.Second()));
      ligand_pocket_orientation_score.First().WriteMDL( outfile);
      io::File::CloseClearFStream( outfile);

      //return score
      return ligand_pocket_orientation_score.Second();
    }

    storage::Pair< FragmentComplete, double> MinimizeLigandInPocket::OrientLigandInPocketEngine
    (
      const FragmentComplete &MOL,
      const FragmentComplete &POCK,
      const size_t &ITERATIONS,
      const size_t &ITERATION_LIMIT,
      const linal::Vector3D &START_POS
    ) const
    {
      FragmentComplete molecule( MOL);
      //BCL_MessageStd( "POS 1");
      //BCL_Debug( molecule.GetCenter());

      //recenter molecule in binding pocket center
      LigandPocketFitScore scorer;
      if( START_POS.IsDefined())
      {
        //BCL_MessageStd( "Scorer option 1");
        molecule.Translate( START_POS - molecule.GetCenter());
        scorer = LigandPocketFitScore( POCK, START_POS);
      }
      else
      {
        //BCL_MessageStd( "Scorer option 2");
        molecule.Translate( POCK.GetCenter() - molecule.GetCenter());
        scorer = LigandPocketFitScore( POCK);
      }

      //BCL_Debug(POCK.GetCenter());
      //BCL_Debug(START_POS);
      //BCL_Debug(molecule.GetCenter());

      //create objects for monte carlo
      MutateMoleculeGeneric mutator( m_Movie);
      util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureAccepted( 0.5, 0.01, ITERATIONS, 1.0, 10));
      mc::Metropolis< double> metropolis( sp_temperature, true);
      opti::CriterionCombine< FragmentComplete, double> criterion_combine;

      // approximator termination criteria
      opti::CriterionNumberIterations< FragmentComplete, double> maximum_number_iterations( ITERATIONS);
      criterion_combine.InsertCriteria(maximum_number_iterations);
      opti::CriterionUnimproved< FragmentComplete, double> unimproved_threshold( ITERATION_LIMIT);
      criterion_combine.InsertCriteria(unimproved_threshold);

      //monte carlo object
      mc::Approximator< FragmentComplete, double> approximator
      (
        scorer,
        mutator,
        metropolis,
        criterion_combine,
        molecule,
        opti::Tracker< FragmentComplete, double>( opti::e_SmallerIsBetter),
        0.0
      );

      // run the approximator
      approximator.Approximate();

      // get the score of the best molecule
      return storage::Pair< FragmentComplete, double>
      (
        approximator.GetTracker().GetBest()->First(),
        approximator.GetTracker().GetBest()->Second()
      );
    }

     //! @brief return parameters for member data that are set up from the labels
     //! @return parameters for member data that are set up from the labels
    io::Serializer MinimizeLigandInPocket::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Orients a small molecule in a pocket cavity by minimizing geometric overlap of matched atoms");
      member_data.AddInitializer
      (
        "iterations",
        "number of iterations for fitting ligand in pocket",
        io::Serialization::GetAgent( &m_Iterations),
        "200"
      );
      member_data.AddInitializer
      (
        "iteration limit",
        "maximum unimproved iterations before monte carlo termination",
        io::Serialization::GetAgent( &m_IterationLimit),
        "40"
      );
      member_data.AddOptionalInitializer
      (
        "reference coordinate",
        "reference coordinate dictating centroid weight",
        io::Serialization::GetAgent( &m_StartPosition)
      );
      return member_data;
    }
    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MinimizeLigandInPocket::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return ConformationComparisonPropertyFieldCorrelation::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace chemistry
} // namespace bcl
