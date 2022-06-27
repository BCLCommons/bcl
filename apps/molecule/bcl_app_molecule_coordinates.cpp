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
#include "bcl_app_molecule_coordinates.h"
#include "chemistry/bcl_chemistry_atom_clash_score.h"
#include "chemistry/bcl_chemistry_fragment_probability_score.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "chemistry/bcl_chemistry_rotamer_dihedral_bond_data.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_isomorphism.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "graph/bcl_graph_undirected_edge.h"
// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    const ApplicationType MoleculeCoordinates::MoleculeCoordinates_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeCoordinates(), GetAppGroups().e_Molecule)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    MoleculeCoordinates::MoleculeCoordinates() :
      m_DihedralBinSizeInDegrees
      (
        new command::FlagStatic
        (
          "dihedral_bin_size",
          "# of degrees each bin will represent in the dihedral angle histograms",
          command::Parameter
          (
            "bin_size",
            "",
            command::ParameterCheckRanged< double>( 0.001, 180.0),
            "30"
          )
        )
      ),
      m_BondAngleBinSizeInDegrees
      (
        new command::FlagStatic
        (
          "bond_angle_bin_size",
          "# of degrees each bin will represent in the bond angle histograms",
          command::Parameter
          (
            "bin_size",
            "",
            command::ParameterCheckRanged< double>( 0.001, 180.0),
            "5"
          )
        )
      ),
      m_StatisticsFlag
      (
        new command::FlagStatic
        (
          "statistics",
          "acquire statistics on bond angles, dihedral angles, and bond lengths"
        )
      ),
      m_DihedralScoresFlag
      (
        new command::FlagStatic
        (
          "dihedral_scores",
          "dihedral scores; "
          "dihedral scores are printed alongside the atom indices of the central bond atoms; "
          "these scores do not include amide bond non-planarity penalties"
        )
      ),
      m_AmideBondDeviationsFlag
      (
        new command::FlagStatic
        (
          "amide_deviations",
          "amide bond deviations from planarity; "
          "amount of deviation from planarity (which is penalized in the BCL::Conf score) "
          "printed alongside the atom indices of the central bond atoms (units: degrees)"
        )
      ),
      m_AmideBondPenaltiesFlag
      (
        new command::FlagStatic
        (
          "amide_penalties",
          "amide bond penalties; "
          "penalty associated with amide bond deviation from planarity; "
          "penalty is applied as the standard weighted sum utilized in BCL::Conf, which is "
          " amide_bond_penalty = amide_deviation_from_planarity( tolerance=10.0 degrees) / 100.0"
          " + amide_deviation_from_planarity( tolerance=15.0 degrees) / 10.0"
          " + amide_deviation_from_planarity( tolerance=25.0 degrees); "
          "printed alongside the atom indices of the central bond atoms (units: ConfScore units)"
        )
      ),
      m_ClashScoresFlag
      (
        new command::FlagStatic
        (
          "clash_scores",
          "atom pair clash scores; "
          "clash scores are printed alongside the atom indices of the clashing atoms; "
          "if hydrogen atoms are to be included in clash score, pass 1; if only heavy "
          "atoms are included in clash score, pass 0 (default)",
          command::Parameter
          (
            "hydrogens",
            "",
            command::ParameterCheckRanged< size_t>( 0, 1),
            "0"
          )
        )
      ),
      m_RecenterMoleculesFlag
      (
        new command::FlagStatic
        (
          "recenter",
          "translate molecules such that the middle of the molecule is at the origin"
        )
      ),
      m_MoleculeCentroidFlag
      (
        new command::FlagStatic
        (
          "centroid",
          "output geometric coordinates for each molecule in ensemble"
        )
      ),
      m_OutputFilenameBase
      (
        new command::FlagStatic
        (
          "output",
          "base name for output of histograms",
          command::Parameter( "output", "base name for output of histograms")
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeCoordinates
    MoleculeCoordinates *MoleculeCoordinates::Clone() const
    {
      return new MoleculeCoordinates( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeCoordinates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeCoordinates::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // add flags for input
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // Dihedral angle size (in degrees)
      sp_cmd->AddFlag( m_DihedralBinSizeInDegrees);

      // Bond angle bin size (in degrees)
      sp_cmd->AddFlag( m_BondAngleBinSizeInDegrees);

      // Statistics
      sp_cmd->AddFlag( m_StatisticsFlag);

      // Dihedral scores
      sp_cmd->AddFlag( m_DihedralScoresFlag);

      // Amide deviations from planarity
      sp_cmd->AddFlag( m_AmideBondDeviationsFlag);

      // Amide nonplanarity penalties
      sp_cmd->AddFlag( m_AmideBondPenaltiesFlag);

      // Atom-pair clash scores
      sp_cmd->AddFlag( m_ClashScoresFlag);

      // Recenter molecules
      sp_cmd->AddFlag( m_RecenterMoleculesFlag);

      // Return centroid coordinates
      sp_cmd->AddFlag( m_MoleculeCentroidFlag);

      // Output filename base
      sp_cmd->AddFlag( m_OutputFilenameBase);

      // molecule reading prefs
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags
      (
        *sp_cmd,
        storage::Set< command::FlagTypeEnum>( command::e_Pthread)
      );

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeCoordinates::GetDescription() const
    {
      return "Performs operations on molecule coordinates";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeCoordinates::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeCoordinates performs operations on coordinates of molecules,\n"
        "Operations include recentering and making histograms of bond lengths, bond and dihedral angles";

      return s_read_me;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeCoordinates::Main() const
    {
      // output streams
      io::OFStream output, output_centroid, output_dihedral_scores, output_amide_deviations, output_amide_penalties, output_clash_scores;

      // output string base
      std::string basename, basename_sdf;
      if( util::EndsWith( io::File::RemoveCompressionExtension( m_OutputFilenameBase->GetFirstParameter()->GetValue()), ".sdf"))
      {
        basename_sdf = m_OutputFilenameBase->GetFirstParameter()->GetValue();
        basename = io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( m_OutputFilenameBase->GetFirstParameter()->GetValue()));
      }
      else
      {
        basename_sdf = m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".sdf";
        basename = m_OutputFilenameBase->GetFirstParameter()->GetValue();
      }

      if( m_RecenterMoleculesFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output,
          basename_sdf
        );
      }
      if( m_MoleculeCentroidFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output_centroid,
          basename + ".centroids.txt"
        );
      }
      if( m_DihedralScoresFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output_dihedral_scores,
          basename + ".dihedral_scores.txt"
        );
      }
      if( m_AmideBondDeviationsFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output_amide_deviations,
          basename + ".amide_deviations.txt"
        );
      }
      if( m_AmideBondPenaltiesFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output_amide_penalties,
          basename + ".amide_penalties.txt"
        );
      }
      if( m_ClashScoresFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output_clash_scores,
          m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".clash_scores.txt"
        );
      }

      // iterate over the fragments from whatever feeds were given
      size_t feed_index( 0);
      if( m_DihedralScoresFlag->GetFlag())
      {
        output_dihedral_scores << "MoleculeIndex,AtomIndexLow,AtomIndexHigh,DihedralScore" << "\n";
      }
      if( m_AmideBondDeviationsFlag->GetFlag())
      {
        output_amide_deviations << "MoleculeIndex,AtomIndexLow,AtomIndexHigh,AmideDeviationFromPlanarity(degrees)" << "\n";
      }
      if( m_AmideBondPenaltiesFlag->GetFlag())
      {
        output_amide_penalties << "MoleculeIndex,AtomIndexLow,AtomIndexHigh,AmideBondPenalty(ConfScore units)" << "\n";
      }
      if( m_ClashScoresFlag->GetFlag())
      {
        output_clash_scores << "MoleculeIndex,AtomIndexLow,AtomIndexHigh,ClashScore" << "\n";
      }
      if( m_MoleculeCentroidFlag->GetFlag())
      {
        output_centroid << "Index,X,Y,Z" << "\n";
      }
      for( chemistry::FragmentFeed feed; feed.NotAtEnd(); ++feed, ++feed_index)
      {
        if( m_StatisticsFlag->GetFlag())
        {
          // add the small molecules' information to all the histograms
          // Skip molecules with bad 3d coordinates
          //if( !feed->HasBadGeometry())
          {
            AddDataToHistograms( *feed);
          }
        }
        if( m_DihedralScoresFlag->GetFlag())
        {
          auto scores( ComputeDihedralScores( *feed));
          WriteScores( output_dihedral_scores, feed_index, scores);
        }
        if( m_AmideBondDeviationsFlag->GetFlag())
        {
          // deviations
          auto scores( feed->GetAmideBondNonPlanarity());
          WriteScores( output_amide_deviations, feed_index, scores);
        }
        if( m_AmideBondPenaltiesFlag->GetFlag())
        {
          // deviations
          auto deviations( feed->GetAmideBondNonPlanarity());

          // initialize tolerances and penalties
          storage::Vector< double> tolerances, penalties;
          tolerances.PushBack( 10.0);
          tolerances.PushBack( 15.0);
          tolerances.PushBack( 25.0);
          penalties.PushBack( 0.01);
          penalties.PushBack( 0.1);
          penalties.PushBack( 1.0);

          // get penalties
          auto amide_penalties( feed->GetPerAmideBondNonPlanarityPenalty( deviations, tolerances, penalties));

          // set deviations to score penalties
          auto scores( deviations);
          for( size_t i( 0); i < deviations.GetSize(); ++i)
          {
            scores( i).Third() = amide_penalties( i);
          }

          WriteScores( output_amide_penalties, feed_index, scores);
        }
        if( m_ClashScoresFlag->GetFlag())
        {
          chemistry::AtomClashScore clash_scorer;
          if( m_ClashScoresFlag->GetFirstParameter()->GetNumericalValue< size_t>() == size_t( 1))
          {
            clash_scorer = chemistry::AtomClashScore( true, false);
          }
          auto scores( clash_scorer.GetClashingPairs( *feed));
          WriteScores( output_clash_scores, feed_index, scores);
        }
        if( m_RecenterMoleculesFlag->GetFlag())
        {
          chemistry::FragmentComplete fragment( *feed);
          fragment.Translate( -fragment.GetCenter());
          fragment.WriteMDL( output);
        }
        if( m_MoleculeCentroidFlag->GetFlag())
        {
          linal::Vector3D centroid( feed->GetCenter());
          WriteCentroid( output_centroid, feed_index, centroid);
        }
      }
      if( m_DihedralScoresFlag->GetFlag())
      {
        io::File::CloseClearFStream( output_dihedral_scores);
      }
      if( m_AmideBondDeviationsFlag->GetFlag())
      {
        io::File::CloseClearFStream( output_amide_deviations);
      }
      if( m_AmideBondPenaltiesFlag->GetFlag())
      {
        io::File::CloseClearFStream( output_amide_penalties);
      }
      if( m_RecenterMoleculesFlag->GetFlag())
      {
        io::File::CloseClearFStream( output);
      }
      if( m_MoleculeCentroidFlag->GetFlag())
      {
        io::File::CloseClearFStream( output_centroid);
      }
      if( m_StatisticsFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output,
          basename + ".dihedral.histograms.txt"
        );
        // Write histogram header
        output << "AtomType1\tBondType12\tAtomType2\tBondType23\tAtomType3\tBondType34\tAtomType4\t";
        WriteHistogram
        (
          output,
          m_DihedralHistogram,
          m_DihedralBinSizeInDegrees->GetFirstParameter()->GetNumericalValue< double>(),
          true
        );
        io::File::CloseClearFStream( output);

        io::File::MustOpenOFStream
        (
          output,
          basename + ".dihedral.sds.txt"
        );
        // Write histogram header
        output << "AtomType1\tBondType12\tAtomType2\tBondType23\tAtomType3\tBondType34\tAtomType4\t";
        WriteStatistics
        (
          output,
          m_BondDihedrals
        );
        io::File::CloseClearFStream( output);

        // open files to write out histograms/stats for bond angles
        io::File::MustOpenOFStream
        (
          output,
          basename + ".bond_angles.stats.txt"
        );
        output << "AtomType1\tBondType12\tAtomType2\tBondType23\tAtomType3\t";
        WriteStatistics( output, m_BondAngleStats);
        io::File::CloseClearFStream( output);

        io::File::MustOpenOFStream
        (
          output,
          basename + ".bond_angles.histograms.txt"
        );
        output << "AtomType1\tBondType12\tAtomType2\tBondType23\tAtomType3\t";
        WriteHistogram
        (
          output,
          m_BondAngleHistogram,
          m_BondAngleBinSizeInDegrees->GetFirstParameter()->GetNumericalValue< double>()
        );
        io::File::CloseClearFStream( output);

        // open files to write out stats for bond lengths
        io::File::MustOpenOFStream
        (
          output,
          basename + ".bond_lengths.stats.txt"
        );
        output << "AtomType1\tBondType12\tAtomType2\t";
        WriteStatistics( output, m_BondLengths);
        io::File::CloseClearFStream( output);
      }
      // end
      return 0;
    }

    //! @brief write histogram given in a map with key identifier to a stream
    //! @param OUTPUT output stream to write to
    //! @param HISTOGRAMS histogram in map with key identifier
    //! @param REFLECTED true if the degrees are reflective (e.g. dihedral bins)
    void MoleculeCoordinates::WriteHistogram
    (
      std::ostream &OUTPUT,
      const storage::Map< std::string, math::Histogram> &HISTOGRAMS,
      const double &BIN_SIZE,
      const bool &REFLECTED
    )
    {
      const size_t number_bins( 180.0 / BIN_SIZE + REFLECTED);
      const bool write_edge_count( BIN_SIZE * int( 180.0 / BIN_SIZE) < 180.0);
      const double offset( REFLECTED ? BIN_SIZE / 2.0 : 0.0);

      // write out bin ranges
      for( size_t bin_number( 0); bin_number < number_bins; ++bin_number)
      {
        OUTPUT << std::max( bin_number * BIN_SIZE - offset, 0.0)
               << '-'
               << std::min( ( bin_number + 1) * BIN_SIZE - offset, 180.0) << '\t';
      }
      if( write_edge_count)
      {
        OUTPUT << BIN_SIZE * number_bins - offset << '-' << 180.0 << '\t';
      }
      OUTPUT << '\n';

      for
      (
        storage::Map< std::string, math::Histogram>::const_iterator
          itr( HISTOGRAMS.Begin()), itr_end( HISTOGRAMS.End());
        itr != itr_end;
        ++itr
      )
      {
        OUTPUT << itr->first << '\t';
        for( size_t bin_number( 0); bin_number < number_bins; ++bin_number)
        {
          OUTPUT << itr->second.GetHistogram()( bin_number) << '\t';
        }
        if( write_edge_count)
        {
          OUTPUT << itr->second.GetBoundariesCounts().Second() << '\t';
        }
        else
        {
          BCL_Assert( itr->second.GetBoundariesCounts().Second() == 0, "Whoops");
        }
        OUTPUT << '\n';
      }
    }

    //! @brief write statistics given in a map with key identifier to a stream
    //! @param OUTPUT output stream to write to
    //! @param STATS statistics in map with key identifier
    void MoleculeCoordinates::WriteStatistics
    (
      std::ostream &OUTPUT,
      const storage::Map< std::string, math::RunningAverageSD< double> > &STATS
    )
    {
      // write out bin description
      OUTPUT << "Count\tMean\tSD\n";

      for
      (
        storage::Map< std::string, math::RunningAverageSD< double> >::const_iterator
          itr( STATS.Begin()), itr_end( STATS.End());
        itr != itr_end;
        ++itr
      )
      {
        OUTPUT << itr->first << '\t'
               << itr->second.GetWeight() << '\t';
        // averages less than 1e-4 are effectively 0; there are not that many significant
        // digits in the coordinates.
        if( itr->second.GetAverage() < 1.0e-4)
        {
          OUTPUT << '0' << '\t';
        }
        else
        {
          OUTPUT << itr->second.GetAverage() << '\t';
        }
        // standard deviations less than 1e-4 are meaningless; there are not that many significant
        // digits in the coordinates.
        if( itr->second.GetStandardDeviation() < 1.0e-4)
        {
          OUTPUT << '0';
        }
        else
        {
          OUTPUT << itr->second.GetStandardDeviation();
        }
        OUTPUT << '\n';
      }
    }

    //! @brief write centroid coordinates to a stream
    //! @param OUTPUT output stream to write to
    //! @param INDEX molecule feed index
    //! @param CENTROID coordinate vector for molecule center
    void MoleculeCoordinates::WriteCentroid
    (
      std::ostream &OUTPUT,
      size_t INDEX,
      const linal::Vector3D CENTROID
    )
    {
      OUTPUT << INDEX << "," << CENTROID.X() << "," << CENTROID.Y() << "," << CENTROID.Z();
      OUTPUT << '\n';
    }

    //! @brief write dihedral scores coordinates to output
    //! @param OUTPUT output stream to write to
    //! @param INDEX molecule feed index
    //! @param SCORES scores and corresponding atom indices
    void MoleculeCoordinates::WriteScores
    (
      std::ostream &OUTPUT,
      size_t INDEX,
      const storage::Vector< storage::Triplet< size_t, size_t, double> > SCORES
    )
    {
      // each score will be its own row; track molecule indices as well so that they can be filtered on column
      for( size_t i( 0), sz( SCORES.GetSize()); i < sz; ++i)
      {
        OUTPUT
        << INDEX << ","
        << SCORES( i).First() << ","
        << SCORES( i).Second() << ","
        << SCORES( i).Third() << '\n';
      }
    }

    //! @brief returns a simple name for a bond type given as a size_t into a string
    //! @returns one of XXXXXX, Aromatic, [1-3]x(Chain|Ring) depending on the bond type
    std::string MoleculeCoordinates::ConvertBondTypeIDToString( const size_t &BOND_TYPE_ID)
    {
      static storage::Vector< std::string> bond_type_names( chemistry::GetConfigurationalBondTypes().GetEnumCount());
      if( BOND_TYPE_ID == chemistry::GetConfigurationalBondTypes().e_Undefined.GetIndex())
      {
        return "XXXXXX";
      }
      else if( bond_type_names( BOND_TYPE_ID).empty())
      {
        chemistry::ConfigurationalBondType bond_type( BOND_TYPE_ID);
        // create the name for this bond type
        if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
        {
          bond_type_names( BOND_TYPE_ID) = "Aromatic";
        }
        else if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Amide)
        {
          bond_type_names( BOND_TYPE_ID) = "Amide";
        }
        else
        {
          // start off with the bond order
          bond_type_names( BOND_TYPE_ID) = util::Format()( bond_type->GetNumberOfElectrons() / 2);
          // followed by ring or chain, as appropriate
          bond_type_names( BOND_TYPE_ID) += bond_type->IsBondInRing() ? "xRing" : "xChain";
        }
      }
      return bond_type_names( BOND_TYPE_ID);
    }

    //! @brief returns a simple name for a bond type given as a size_t into a string
    //! @returns one of XXXXXX, Aromatic, Single, Double, or Triple depending on the bond type
    std::string MoleculeCoordinates::ConvertBondTypeIDToStringIgnoreTopology( const size_t &BOND_TYPE_ID)
    {
      static storage::Vector< std::string> bond_type_names( chemistry::GetConfigurationalBondTypes().GetEnumCount());
      if( BOND_TYPE_ID == chemistry::GetConfigurationalBondTypes().e_Undefined.GetIndex())
      {
        return "XXXXXX";
      }
      else if( bond_type_names( BOND_TYPE_ID).empty())
      {
        chemistry::ConfigurationalBondType bond_type( BOND_TYPE_ID);
        // create the name for this bond type
        if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Aromatic)
        {
          bond_type_names( BOND_TYPE_ID) = "Aromatic";
        }
        else if( bond_type->GetConjugation() == chemistry::ConstitutionalBondTypeData::e_Amide)
        {
          bond_type_names( BOND_TYPE_ID) = "Amide";
        }
        else if( bond_type->GetNumberOfElectrons() == size_t( 2))
        {
          bond_type_names( BOND_TYPE_ID) = "Single";
        }
        else if( bond_type->GetNumberOfElectrons() == size_t( 4))
        {
          bond_type_names( BOND_TYPE_ID) = "Double";
        }
        else if( bond_type->GetNumberOfElectrons() == size_t( 6))
        {
          bond_type_names( BOND_TYPE_ID) = "Triple";
        }
        else
        {
          bond_type_names( BOND_TYPE_ID) = "XXXXXX";
        }
      }
      return bond_type_names( BOND_TYPE_ID);
    }

    //! @brief add the data from a small molecule to the histograms
    //! @param MOLECULE the small molecule whose data to add
    void MoleculeCoordinates::AddDataToHistograms( const chemistry::ConformationInterface &MOLECULE) const
    {
      static const double dihedral_bin_size
      (
        m_DihedralBinSizeInDegrees->GetFirstParameter()->GetNumericalValue< double>()
      );
      static const double bondangle_bin_size
      (
        m_BondAngleBinSizeInDegrees->GetFirstParameter()->GetNumericalValue< double>()
      );
      static const math::Histogram empty_dihedral_histogram( -dihedral_bin_size / 2.0, dihedral_bin_size, 1 + 180.0 / dihedral_bin_size);
      static const math::Histogram empty_bondangle_histogram( 0.0, bondangle_bin_size, 180.0 / bondangle_bin_size);

      storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> >
        dihedrals( MOLECULE.GetDihedralAnglesByType());

      // add the small molecule's dihedral angle info to the histograms
      for
      (
        storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> >::const_iterator
          itr( dihedrals.Begin()), itr_end( dihedrals.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !ContainsOnlyGasteigerAtomTypes( itr->first))
        {
          continue;
        }
        std::string angle_string( ConvertAtomBondTypeIntoString( itr->first));
        if( !m_DihedralHistogram.Has( angle_string))
        {
          m_DihedralHistogram[ angle_string] = empty_dihedral_histogram;
        }
        math::Histogram &histogram( m_DihedralHistogram[ angle_string]);
        auto &avesd( m_BondDihedrals[ angle_string]);
        for
        (
          storage::Vector< double>::const_iterator
            itr_angles( itr->second.Begin()), itr_angles_end( itr->second.End());
          itr_angles != itr_angles_end;
          ++itr_angles
        )
        {
          double degree( math::Angle::Degree( *itr_angles));
          histogram.PushBack( degree);
          if( degree < 15.0)
          {
            degree += 360.0;
          }
          // bins are centered about 0, 30, etc. So standard deviation of the bins relative to the bin centers
          avesd += degree - 30.0 * int( ( degree - 15.0) / 30.0 + 1.0);
        }
      }

      storage::Map< storage::VectorND< 5, size_t>, storage::Vector< double> >
        angles( MOLECULE.GetBondAnglesByType());

      // add the small molecule's dihedral angle info to the histograms
      for
      (
        storage::Map< storage::VectorND< 5, size_t>, storage::Vector< double> >::const_iterator
          itr( angles.Begin()), itr_end( angles.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !ContainsOnlyGasteigerAtomTypes( itr->first))
        {
          continue;
        }
        std::string angle_string( ConvertAtomBondTypeIntoString( itr->first));
        if( !m_BondAngleHistogram.Has( angle_string))
        {
          m_BondAngleHistogram[ angle_string] = empty_bondangle_histogram;
        }
        math::Histogram &histogram( m_BondAngleHistogram[ angle_string]);
        math::RunningAverageSD< double> &statistic( m_BondAngleStats[ angle_string]);
        for
        (
          storage::Vector< double>::const_iterator
            itr_angles( itr->second.Begin()), itr_angles_end( itr->second.End());
          itr_angles != itr_angles_end;
          ++itr_angles
        )
        {
          histogram.PushBack( math::Angle::Degree( *itr_angles));
          statistic += math::Angle::Degree( *itr_angles);
        }
      }

      storage::Map< storage::VectorND< 3, size_t>, storage::Vector< double> >
        lengths( MOLECULE.GetBondLengthsByType());
      // Ring and chain bonds have identical statistics, except aromatic bonds, so replace accordingly

      // add the small molecule's lengths to the histograms
      for
      (
        storage::Map< storage::VectorND< 3, size_t>, storage::Vector< double> >::const_iterator
          itr( lengths.Begin()), itr_end( lengths.End());
        itr != itr_end;
        ++itr
      )
      {
        // create atom types out of the information
        chemistry::AtomType atom_type_a( itr->first.First());
        chemistry::AtomType atom_type_b( itr->first.Third());
        chemistry::AtomType base_atom_type_a( atom_type_a->GetBaseAtomType());
        chemistry::AtomType base_atom_type_b( atom_type_b->GetBaseAtomType());
        std::string base_atom_type_a_name( base_atom_type_a->GetElementType()->GetChemicalSymbol());
        if( base_atom_type_a->GetFormalCharge() >= 0)
        {
          base_atom_type_a_name += '+';
        }
        base_atom_type_a_name += util::Format()( base_atom_type_a->GetFormalCharge());
        std::string base_atom_type_b_name( base_atom_type_b->GetElementType()->GetChemicalSymbol());
        if( base_atom_type_b->GetFormalCharge() >= 0)
        {
          base_atom_type_b_name += '+';
        }
        base_atom_type_b_name += util::Format()( base_atom_type_b->GetFormalCharge());

        // get the basic bond type string
        const std::string bond_type( ConvertBondTypeIDToStringIgnoreTopology( itr->first.Second()));

        util::SiPtrVector< math::RunningAverageSD< double> > bond_length_stats;
        // create a siptr vector of statistics to add these bond statistics to
        // bond stats will be added at up to 4 levels:
        // AtomType1 BondType AtomType2
        // AtomType1 BondType ElementTypeCharge2
        // AtomType2 BondType ElementTypeCharge1
        // ElementTypeCharge1 BondType ElementTypeCharge2
        // get the base string
        std::string detailed_type( base_atom_type_a_name + ' ' + bond_type + ' ' + base_atom_type_b_name);
        bond_length_stats.PushBack( &m_BondLengths[ detailed_type]);
        if( atom_type_a->IsGasteigerAtomType())
        {
          detailed_type = atom_type_a.GetName() + ' ' + bond_type + ' ' + base_atom_type_b_name;
          bond_length_stats.PushBack( &m_BondLengths[ detailed_type]);
          if( atom_type_b->IsGasteigerAtomType())
          {
            detailed_type = atom_type_a.GetName() + ' ' + bond_type + ' ' + atom_type_b.GetName();
            bond_length_stats.PushBack( &m_BondLengths[ detailed_type]);
          }
        }
        if( atom_type_b->IsGasteigerAtomType())
        {
          detailed_type = atom_type_b.GetName() + ' ' + bond_type + ' ' + base_atom_type_a_name;
          bond_length_stats.PushBack( &m_BondLengths[ detailed_type]);
        }

        // add bond length stats
        for
        (
          util::SiPtrVector< math::RunningAverageSD< double> >::iterator
            itr_stats( bond_length_stats.Begin()), itr_stats_end( bond_length_stats.End());
          itr_stats != itr_stats_end;
          ++itr_stats
        )
        {
          math::RunningAverageSD< double> &statistic( **itr_stats);
          for
          (
            storage::Vector< double>::const_iterator
              itr_angles( itr->second.Begin()), itr_angles_end( itr->second.End());
            itr_angles != itr_angles_end;
            ++itr_angles
          )
          {
            statistic += *itr_angles;
          }
        }
      }
    }

    //! @brief compute the dihedral angle scores for the molecule
    //! @param MOLECULE the small molecule whose dihedral scores will be computed
    //! @return a vector of scores and the atom indices of the central dihedral bond atoms
    storage::Vector< storage::Triplet< size_t, size_t, double> > MoleculeCoordinates::ComputeDihedralScores
    (
      const chemistry::ConformationInterface &MOLECULE
//      const chemistry::ConformationInterface &molecule
    )
    {

      // note 09-09-2021
      // FYI for direct comparison to bcl.exe molecule:ConformerGenerator -max_iterations 0, requires standardizing bond lengths and updating hydrogens
//      chemistry::FragmentComplete MOLECULE( molecule);
//      MOLECULE.StandardizeBondLengths();
//      MOLECULE.UpdateH();

      // search fragments for molecule of interest
      static chemistry::RotamerLibraryFile rotamer_library;
      static util::ShPtr< chemistry::SearchFragmentLibraryFromTree> search_fragment_library( new chemistry::SearchFragmentLibraryFromTree( rotamer_library));
      util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> fragments_iso_lib( search_fragment_library->FindFragmentsOfMolecule( MOLECULE));

      util::ShPtrVector< chemistry::RotamerDihedralBondData> bond_mapping
      (
        chemistry::SmallMoleculeFragmentMapping().MapFragmentIsomorphisms( MOLECULE, fragments_iso_lib)
      );

      // This is how you get the dihedral bond central bond atom pair and bond type
      storage::Vector< graph::UndirectedEdge< chemistry::ConfigurationalBondType> > dihedral_edges;
      if( !bond_mapping.IsEmpty())
      {
        dihedral_edges = chemistry::PriorityDihedralAngles::GetDihedralEdges( bond_mapping.FirstElement()->GetFragment().GetMolecule());
      }

      // compute the unweighted dihedral scores
      chemistry::FragmentProbabilityScore scorer( bond_mapping, false); // should probably be a flag, though not sure the statistics in stereoisomers warrant it
      storage::Vector< double> dihedral_scores( scorer.GetDihedralScores( MOLECULE));
      storage::Vector< size_t> dihedral_component_sizes( scorer.GetDihedralComponentSize());
      size_t number_dihedral_bonds( scorer.GetNumberDihedralBonds());

      // normalize dihedral scores by component sizes
      double total_weight( 0.0);
      for( size_t i( 0); i < number_dihedral_bonds; ++i)
      {
        if( dihedral_scores( i))
        {
          total_weight += 1.0 / double( dihedral_component_sizes( i));
        }
      }

      for( size_t i( 0); i < number_dihedral_bonds; ++i)
      {
        // checks for 0-weighted positions to ignore them
        if( dihedral_scores( i))
        {
          dihedral_scores( i) *= 1.0 / double( dihedral_component_sizes( i)) / total_weight;
        }
      }

      // these should be the same size
      BCL_Assert
      (
        dihedral_edges.GetSize() == dihedral_scores.GetSize(),
        "Number of dihedral edges does not match the number of dihedral scores!!"
      );

      // get the edge indices (atoms contributing to the edge)
      storage::Vector< storage::Triplet< size_t, size_t, double> > output( dihedral_scores.GetSize());
      size_t score_index( 0);
      for
      (
          auto edge_itr( dihedral_edges.Begin()), edge_itr_end( dihedral_edges.End());
          edge_itr != edge_itr_end;
          ++edge_itr, ++score_index
      )
      {
        output( score_index) =
            storage::Triplet< size_t, size_t, double>
        (
          edge_itr->GetIndexLow(),
          edge_itr->GetIndexHigh(),
          dihedral_scores( score_index)
        );
      }
      return output;
    }

  } // namespace app
} // namespace bcl

