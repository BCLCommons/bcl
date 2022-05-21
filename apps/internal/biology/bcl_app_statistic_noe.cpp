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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_statistics.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "score/bcl_score_restraint_noe_knowledge_based.h"

namespace bcl
{
  namespace app
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StatisticNoe
    //! @brief app is for creating NOE distance statistics between HH and Cb. This includes HA and HN and all are traced
    //!        back to the Cb
    //!
    //! @author akinlr, weinerbe
    //! @date 01/03/11
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StatisticNoe :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! list of pdbs for doing statistics over
      util::ShPtr< command::FlagInterface> m_PDBListFlag;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_HistogramOutputFileFlag;

      //! flag for specifying maximum NOE distance length
      util::ShPtr< command::FlagInterface> m_DistanceCutoffFlag;

      //! flag for specifying minimum number of residues between AA
      util::ShPtr< command::FlagInterface> m_ResidueDistanceFlag;

      //! flag for specifying whether or not to write histograms
      util::ShPtr< command::FlagInterface> m_WriteHistogramsFlag;

      //! flag for specifying whether to use H and HA atoms
      util::ShPtr< command::FlagInterface> m_BackboneProtonFlag;

      //! flag for writing out individual potentials
      util::ShPtr< command::FlagInterface> m_RosettaPotentials;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StatisticNoe();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      StatisticNoe *Clone() const
      {
        return new StatisticNoe( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief GetHydrogenAtoms gives a list of Side Chain atoms
      //! @param PROTEIN_MODEL the list from which the hydrogen atoms will come
      //! @return returns a List of hydrogen atoms paired with their AA
      storage::List
      <
        storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
      > GetHydrogenAtoms( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief GetHydrogenAtomPairs gives a list of pairs of hydrogen atoms within a single protein
      //! @param HYDROGEN_ATOM_LIST the list of all possible hydrogen atoms which can be paired
      //! @return returns a List of pairs of hydrogen atoms and their AA that are at least 8A apart and are no closer
      //!         than 6 amino acids in sequence
      storage::List
      <
        storage::Pair
        <
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >,
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
        >
      > GetHydrogenAtomPairs
      (
        const storage::List
        <
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
        > &HYDROGEN_ATOM_LIST
      ) const;

      //! @brief GetProteinModel creates a protein model from a pdb filename
      //! @param PDB_FILENAME is the name of the pdb filename from which a protein model will be created
      //! @return assemble::ProteinModel created from "PDB_FILENAME"
      assemble::ProteinModel GetProteinModel( const std::string &PDB_FILENAME) const;

      //! @brief GetNumberOfBonds determines how many bonds are between the hydrogen of a side chain atom and the Cb
      //! @param ATOM_STRING string representing the atom type
      //! @return number of bonds between Cb and hydrogen atom
      size_t GetNumberOfBonds( const std::string &ATOM_STRING) const;

      //! @brief returns whether a hydrogen atom type is in the backbone
      //! @return whether a hydrogen atom type is in the backbone
      bool IsBackboneProton( const biol::AtomType &ATOM_TYPE) const;

      //! @brief returns whether a hydrogen atom type is H
      //! @return whether a hydrogen atom type is H
      bool IsBackboneH( const biol::AtomType &ATOM_TYPE) const;

      //! @brief returns whether a hydrogen atom type is HA
      //! @return whether a hydrogen atom type is HA
      bool IsBackboneHA( const biol::AtomType &ATOM_TYPE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType StatisticNoe_Instance;

    }; // StatisticNoe

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> StatisticNoe::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flag to input list of pdbs
      sp_cmd->AddFlag( m_PDBListFlag);

      // add flag to input where to put the histogram output
      sp_cmd->AddFlag( m_HistogramOutputFileFlag);

      // add flag for distance cutoff
      sp_cmd->AddFlag( m_DistanceCutoffFlag);

      // add flag for sequence distance minimum
      sp_cmd->AddFlag( m_ResidueDistanceFlag);

      // add flag for writing histograms
      sp_cmd->AddFlag( m_WriteHistogramsFlag);

      // add flag for using backbone protons
      sp_cmd->AddFlag( m_BackboneProtonFlag);

      sp_cmd->AddFlag( m_RosettaPotentials);

      // pdb factory flags
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());
      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter( biol::GetAAClasses().e_AAComplete);
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int StatisticNoe::Main() const
    {
      // if writing out current histogram as individual potentials
      if( m_RosettaPotentials->GetFlag())
      {
        // create IFStream "read"
        io::IFStream read;

        // open "read" and bind it to the histogram file containing SL-CB distances
        io::File::MustOpenIFStream
        (
          read,
          score::Score::AddHistogramPath( score::RestraintNoeKnowledgeBased::GetDefaultHistogramFilename())
        );

        // initialize atom type
        biol::AtomType atom_type;

        // iterate through the entire histogram
        while( read >> atom_type && !read.eof())
        {
          // break if undefined atom_type (end of file)
          if( atom_type == biol::GetAtomTypes().e_Undefined)
          {
            break;
          }

          // create an empty size_t for storing the bond distances
          size_t bond_distance;

          // read histogram
          math::Histogram histogram;
          read >> bond_distance >> histogram;

          if( bond_distance > 0)
          {
            // create math::Vector values and initialize with the values in the histogram provided
            linal::Vector< double> values( histogram.GetHistogram());

            // add pseudocount
            const double pseudocount( 100);
            values += pseudocount;

            // variable that adds up all counts
            const double total_sum( values.Sum());

            // convert to propensities
            values *= values.GetSize() / total_sum;

            //-log of every bin to get the energy
            for( double *ptr( values.Begin()), *ptr_end( values.End()); ptr != ptr_end; ++ptr)
            {
              *ptr = -log( *ptr);
            }

            // need to normalize the energies of "values" to be between 0 and 1
            values /=
              (
                math::Statistics::MaximumValue( values.Begin(), values.End()) -
                math::Statistics::MinimumValue( values.Begin(), values.End())
              );

            // need to shift the energy potential to be 0 or less (i.e. no penalties)
            values -= math::Statistics::MaximumValue( values.Begin(), values.End());

            // write out the histogram
            const std::string filename( atom_type.GetName() + "_" + util::Format()( bond_distance) + ".potential");
            io::OFStream write;
            io::File::MustOpenOFStream( write, filename);

            // write out the bins
            write << "x_axis\t";
            const size_t nr_bins( histogram.GetNumberOfBins());
            for( size_t i( 0); i != nr_bins; ++i)
            {
              const double bin_size( histogram.GetBinSize());
              const double bin( histogram.GetBoundaries().First() + double( i) * bin_size + 0.5 * bin_size);
              write << bin << '\t';
            }

            // write out the values
            write << "\ny_axis\t";
            for( double *ptr( values.Begin()), *ptr_end( values.End()); ptr != ptr_end; ++ptr)
            {
              write << *ptr << '\t';
            }

            io::File::CloseClearFStream( write);
          }
        }

        // close and clear read stream
        io::File::CloseClearFStream( read);
        return 0;
      }

      // create Vector "pdb_list" initialize with the list of strings obtained from
      // contained in the file "m_PDBListFlag"
      const storage::Vector< std::string> pdb_list
      (
        command::StringVectorFromFilenameParameter( *( m_PDBListFlag->GetFirstParameter()))
      );

      // create a vector to store data for histogram of side chain distance + bond # -Cb distance
      storage::Vector< double> histogram_knowledge_based_data;

      // create a map to store data for 2D histogram if needed
      storage::Map< size_t, storage::Vector< storage::VectorND< 2, double> > > histogram2d_data;

      // create a map to store data for 1D histogram
      storage::Map
      <
        storage::Pair< biol::AtomType, size_t>, storage::Vector< double>
      > histogram_data;

      // loop through "file_path_and_name" to do statistics over all proteins
      for
      (
        storage::Vector< std::string>::const_iterator file_itr( pdb_list.Begin()), file_itr_end( pdb_list.End());
        file_itr != file_itr_end;
        ++file_itr
      )
      {
        // create ProteinModel "protein_model" and initialize with protein model created from "pdb_filename"
        const assemble::ProteinModel protein_model( GetProteinModel( *file_itr));

        // message that the protein model "pdb_filename" was created
        BCL_MessageCrt( "protein model " + *file_itr + " created");

        // create a list of hydrogen side chain atoms in the protein for each model
        const storage::List
        <
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
        > pdb_hyd_list( GetHydrogenAtoms( protein_model));

        // use the list of hydrogen side chain atoms to create pairs of these atoms
        const storage::List
        <
          storage::Pair
          <
            storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >,
            storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
          >
        > pdb_hyd_pair_list( GetHydrogenAtomPairs( pdb_hyd_list));

        // iterate over the hydrogen atom pairs
        for
        (
          storage::List
          <
            storage::Pair
            <
              storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >,
              storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
            >
          >::const_iterator pair_itr( pdb_hyd_pair_list.Begin()), pair_itr_end( pdb_hyd_pair_list.End());
          pair_itr != pair_itr_end; ++pair_itr
        )
        {
          // calculate the distance between the hydrogen atoms
          const double hyd_atom_distance
          (
            linal::Distance( pair_itr->First().First()->GetCoordinates(), pair_itr->Second().First()->GetCoordinates())
          );

          // get the coordinates for the main heavy atoms
          linal::Vector3D coords_a( pair_itr->First().Second()->GetFirstSidechainAtom().GetCoordinates());
          linal::Vector3D coords_b( pair_itr->Second().Second()->GetFirstSidechainAtom().GetCoordinates());

          // if using backbone protons
          if( m_BackboneProtonFlag->GetFlag())
          {
            // if the first atom type is a backbone proton
            if( IsBackboneProton( pair_itr->First().First()->GetType()))
            {
              // set the coords
              coords_a = pair_itr->First().First()->GetCoordinates();
            }
            // if the second atom type is a backbone proton
            if( IsBackboneProton( pair_itr->Second().First()->GetType()))
            {
              // set the coords
              coords_b = pair_itr->Second().First()->GetCoordinates();
            }
          }

          // calculate the distance between the heavy atoms of choice
          const double main_atom_distance( linal::Distance( coords_a, coords_b));

          // calculate how many bonds are between the hydrogen atoms and their Cb
          const size_t nr_bonds_a( GetNumberOfBonds( pair_itr->First().First()->GetType().GetName()));
          const size_t nr_bonds_b( GetNumberOfBonds( pair_itr->Second().First()->GetType().GetName()));

          // if either bond distance is undefined
          if( !util::IsDefined( nr_bonds_a) || !util::IsDefined( nr_bonds_b))
          {
            // skip this pair
            continue;
          }

          const size_t bond_difference
          (
            math::Absolute
            (
              GetNumberOfBonds( pair_itr->First().First()->GetType().GetName()) +
              GetNumberOfBonds( pair_itr->Second().First()->GetType().GetName())
            )
          );

          if( m_WriteHistogramsFlag->GetFlag())
          {
            // create a vector containing the two distances calculated
            histogram2d_data[ bond_difference].PushBack
            (
              storage::VectorND< 2, double>( main_atom_distance, hyd_atom_distance)
            );

            // create a vector containing for 1D histogram of Hyd distance + bond number - Cb distance
            histogram_knowledge_based_data.PushBack( hyd_atom_distance + bond_difference - main_atom_distance);
          }

          // set the atom type to CB
          biol::AtomType atom_type( biol::GetAtomTypes().CB);

          // if backbone protons are to be used
          if( m_BackboneProtonFlag->GetFlag())
          {
            const biol::AtomType &atom_a_type( pair_itr->First().First()->GetType());
            const biol::AtomType &atom_b_type( pair_itr->Second().First()->GetType());

            // if either atom type is H
            if( IsBackboneH( atom_a_type) || IsBackboneH( atom_b_type))
            {
              atom_type = biol::GetAtomTypes().H;
            }
            // if either atom type is HA
            if( IsBackboneHA( atom_a_type) || IsBackboneHA( atom_b_type))
            {
              atom_type = biol::GetAtomTypes().HA;
            }
          }

          // create a vector for the 1D histogram
          histogram_data
          [
            storage::Pair< biol::AtomType, size_t>( atom_type, bond_difference)
          ].PushBack( hyd_atom_distance - main_atom_distance);

        } // iteration through hydrogen atom pairs
       } // iteration through each protein

      // create a 1D histogram of SC distance - Cb distance
      math::Histogram histogram( -15.0, 0.25, 80);

      // write "histogram" to file designated by "m_HistogramOutputFilename"
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_HistogramOutputFileFlag->GetFirstParameter()->GetValue());

      // iterate through the histogram map to create 1D histograms for each bond_difference
      for
      (
        storage::Map
        <
          storage::Pair< biol::AtomType, size_t>,
          storage::Vector< double>
        >::const_iterator map_itr( histogram_data.Begin()), map_itr_end( histogram_data.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        // create a 1D histogram of SC distance - Cb distance
        math::Histogram this_histogram( -15.0, 0.25, 80);

        // calculate the 1D histogram
        this_histogram.CalculateHistogram( map_itr->second);
        this_histogram.RemoveBinsAfterIndex
        (
          std::min( this_histogram.GetIndexOfLastInformationContainingBin() + 1, this_histogram.GetNumberOfBins())
        );
        this_histogram.RemoveBinsBeforeIndex
        (
          std::max( 0, int( this_histogram.GetIndexOfFirstInformationContainingBin()) - 1)
        );

        write << map_itr->first.First() << '\n';
        write << map_itr->first.Second() << '\n';
        write << this_histogram << '\n';
      }

      io::File::CloseClearFStream( write);

      if( m_WriteHistogramsFlag->GetFlag())
      {
        // create and calculate a 1D histogram for a SC distance + bond numbers - Cb distance
        math::Histogram histogram_knowledge_based( -30.0, 1.0, 50);
        histogram_knowledge_based.CalculateHistogram( histogram_knowledge_based_data);

        // create an ofstream for linear gnuplot file
        io::OFStream write;
        std::string title_final( "histogram_knowledge_based");
        io::File::MustOpenOFStream( write, title_final + ".plot");

        //write the plot
        histogram_knowledge_based.WriteLinearGnuplot( write, title_final);

        // close the stream
        io::File::CloseClearFStream( write);

        for
        (
            storage::Map
            <
              storage::Pair< biol::AtomType, size_t>,
              storage::Vector< double>
            >::const_iterator map_itr( histogram_data.Begin()), map_itr_end( histogram_data.End());
          map_itr != map_itr_end;
          ++map_itr
        )
        {
          // create an ofstream for gnuplot file
          io::OFStream write;
          std::string title( "histogram_" + map_itr->first.First().GetName() + util::Format()( map_itr->first.Second()));
          io::File::MustOpenOFStream( write, title + ".plot");

          histogram.CalculateHistogram( map_itr->second);

          // write the plot
          histogram.WriteLinearGnuplot( write, title);

          // close the stream
          io::File::CloseClearFStream( write);
        }

        // create a 2D histogram
        math::Histogram2D histogram_2d
        (
          storage::VectorND< 2, double>( 0.0, 0.0),
          storage::VectorND< 2, double>( 1.0, 1.0),
          storage::VectorND< 2, size_t>( 16, 8)
        );

        // iterate through the histogram map to create 2D histograms for each bond_difference
        for
        (
          storage::Map
          <
            size_t, storage::Vector< storage::VectorND< 2, double> >
          >::const_iterator map_itr( histogram2d_data.Begin()), map_itr_end( histogram2d_data.End());
          map_itr != map_itr_end;
          ++map_itr
        )
        {
          // calculate the 2D histogram
          histogram_2d.CalculateHistogram( map_itr->second);

          // create an ofstream for gnuplot file
          io::OFStream write;
          std::string title_2d( "histogram2D_" + util::Format()( map_itr->first));
          io::File::MustOpenOFStream( write, title_2d + ".plot");

          // write the plot
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromHistogram( histogram_2d, true, true);
          heatmap.SetTitleAndLabel( title_2d, "cb-distance", "noe", "count");
          heatmap.WriteScript( write);

          // close the stream
          io::File::CloseClearFStream( write);
        } // end of creating 2D histogram
      } // end of creating histograms if flag is set

      // end
      return 0;
     }

    // default constructor
    StatisticNoe::StatisticNoe() :
      m_PDBListFlag
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\tlist of pdbs to do statistics over",
          command::Parameter
          (
            "string",
            "\tname of file which is the list of pdbs to do statistics over",
            ""
          )
        )
      ),
      m_HistogramOutputFileFlag
      (
        new command::FlagStatic
        (
          "output_file",
          "\tpath and name of the output file which will hold the results.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "noe_knowledge_based.histogram"
          )
        )
      ),
      m_DistanceCutoffFlag
      (
        new command::FlagStatic
        (
          "distance_cutoff",
          "\tthe maximum allowed distance considered",
          command::Parameter
          (
            "distance_value",
            "\tthe maximum allowed distance considered",
            "6.0"
          )
        )
      ),
      m_ResidueDistanceFlag
      (
        new command::FlagStatic
        (
          "residue_distance",
          "\tthe minimum sequence distance between two residues",
          command::Parameter
          (
            "distance_value",
            "\tthe minimum sequence distance between two residues",
            "5"
          )
        )
      ),
      m_WriteHistogramsFlag
      (
        new command::FlagStatic
        (
          "write_histogram",
          "\twrite several histograms pertaining to different distance statistics of protein list",
          command::Parameter
          (
            "write_histogram",
            "\twrite several histograms pertaining to distance statistics of protein list",
            ""
          )
        )
      ),
      m_BackboneProtonFlag
      (
        new command::FlagStatic
        (
          "backbone_protons",
          "\tcalculate distances to backbone protons rather than using the CB coordinates for those "
            "atom types (H and HA)"
        )
      ),
      m_RosettaPotentials
      (
        new command::FlagStatic
        (
          "rosetta_potentials",
          "\twrite individual potentials calculated from the current histogram file to be used by Rosetta"
        )
      )
    {
    }

    const ApplicationType StatisticNoe::StatisticNoe_Instance
    (
      GetAppGroups().AddAppToGroup( new StatisticNoe(), GetAppGroups().e_Restraint)
    );

    //! @brief GetHydrogenAtoms gives a list of Hydrogen atoms from a protein model
    //! @param PROTEIN_MODEL the list of all possible amino acids whose side chain atoms can be used
    //! @return returns a list of hydrogen atoms and the amino acid they came from
    storage::List
    <
      storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
    > StatisticNoe::GetHydrogenAtoms
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create an empty list to store hydrogen atom and AA data in
      storage::List
      <
        storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
      > hydrogen_atom_list;

      // create a vector of the amino acids from the protein model
      util::SiPtrVector< const biol::AABase> aa_vector( PROTEIN_MODEL.GetAminoAcids());

      // iterate through the vector of amino acids
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( aa_vector.Begin()),
          aa_itr_end( aa_vector.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // iterate through the atoms in the amino acids
        for
        (
          util::SiPtrVector< const biol::Atom>::const_iterator atom_itr( ( *aa_itr)->GetAtoms().Begin()),
            atom_itr_end( ( *aa_itr)->GetAtoms().End());
          atom_itr != atom_itr_end;
          ++atom_itr
        )
        {
          // store only the hydrogen atoms that are on the side chain in a pair with the original amino acid
          if( ( *atom_itr)->GetType()->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
          {
            // store atoms that qualify
            hydrogen_atom_list.PushBack
            (
              storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >( *atom_itr, *aa_itr)
            );
          }
        } // end iteration through atoms in aas
      } // end iteration through vector of aa

      return hydrogen_atom_list;
    }

    //! @brief GetHydrogenAtomPairs gives a list of pairs of hydrogen atoms
    //! @param HYDROGEN_ATOM_LIST the list of all possible hydrogen atoms to pair
    //! @return returns a List with pairs of hydrogen atoms and their respective AA
    storage::List
    <
      storage::Pair
      <
        storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >,
        storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
      >
    > StatisticNoe::GetHydrogenAtomPairs
    (
      const storage::List
      <
        storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
      > &HYDROGEN_ATOM_LIST
    ) const
    {
      // create an empty list to contain a pair of atoms that are paired with their respective AA
      storage::List
      <
        storage::Pair
        <
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >,
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
        >
      > hydrogen_atom_pairs;

      // iterate through the list of hydrogen atoms
      for
      (
        storage::List
        <
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
        >::const_iterator hyd_atm_itr_a( HYDROGEN_ATOM_LIST.Begin()), hyd_atm_itr_end( HYDROGEN_ATOM_LIST.End());
        hyd_atm_itr_a != hyd_atm_itr_end;
        ++hyd_atm_itr_a
      )
      {
        // iterate through the same list a second time in order to create a new list of paired atoms
        storage::List
        <
          storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
        >::const_iterator hyd_atm_itr_b( hyd_atm_itr_a);
        ++hyd_atm_itr_b;

        for( ; hyd_atm_itr_b != hyd_atm_itr_end; ++hyd_atm_itr_b)
        {
          // include atoms as pairs only if they are equal to or further apart than m_DistanceCutoffFlag
          // and only include atoms as pairs if they are as close as the m_ResidueDistanceFlag in space
          if
          (
            size_t( math::Absolute( hyd_atm_itr_a->Second()->GetSeqID() - hyd_atm_itr_b->Second()->GetSeqID())) >=
            m_DistanceCutoffFlag->GetFirstParameter()->GetNumericalValue< size_t>() &&
            linal::Distance( hyd_atm_itr_a->First()->GetCoordinates(), hyd_atm_itr_b->First()->GetCoordinates()) <=
            m_ResidueDistanceFlag->GetFirstParameter()->GetNumericalValue< double>()
          )
          {
            // store the atom pairs in a new list
            hydrogen_atom_pairs.PushBack
            (
              storage::Pair
              <
                storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >,
                storage::Pair< util::SiPtr< const biol::Atom>, util::SiPtr< const biol::AABase> >
              >
              (
                *hyd_atm_itr_a, *hyd_atm_itr_b
              )
            );
          } // if statement pair storage
        } // second iteration
      } // first iteration

      return hydrogen_atom_pairs;
    } // end GetHydrogenAtomPairs

    //! @brief GetProteinModel is for creating a protein model given a pdb filename
    //! @param PDB_FILENAME is the name of the pdb filename from which a protein model will be created
    //! @return assemble::ProteinModel created from "PDB_FILENAME"
    assemble::ProteinModel StatisticNoe::GetProteinModel( const std::string &PDB_FILENAME) const
    {
      // output of current pdb filename
      BCL_MessageStd( "processing pdb: " + PDB_FILENAME);

      // create static PDBFactory "pdb_factory" and initialize it to work with AABACKBONE type amino acids
      static const pdb::Factory pdb_factory;

      // create io::IFStream "read"
      io::IFStream read;

      // open "read" and bind it to "pdb_filename"
      io::File::MustOpenIFStream( read, PDB_FILENAME);

      // create PDBReader "pdb" from "read"
      pdb::Handler pdb( read);

      // close and clear "read"
      io::File::CloseClearFStream( read);

      // message telling which pdb is being processed
      BCL_MessageStd( "processing pdb: " + PDB_FILENAME + " is complete");

      // create ProteinModel "protein_model" and initialize with "pdb"
      const assemble::ProteinModel protein_model( pdb_factory.ProteinModelFromPDB( pdb));

      // return "protein_model"
      return protein_model;
    }

    //! @brief determines how many bonds are between a hydrogen atom on a side chain and the Cb
    //! @param ATOM_STRING string representing the atom type
    //! @return number of bonds between Cb and sidechain atom's hydrogen
    size_t StatisticNoe::GetNumberOfBonds( const std::string &ATOM_STRING) const
    {

      // read in the second character to determine the side chain position (i.e. beta, gamma, delta, etc.)
      const char atom_pos( ATOM_STRING[ 1]);

      // initialize number of bonds between proton and CB
      size_t nr_bonds( 0);

      // if backbone protons are to be used directly and the type is right
      if
      (
        m_BackboneProtonFlag->GetFlag() &&
        ( ATOM_STRING == "H" || ATOM_STRING == "HA" || ATOM_STRING == "HA2" || ATOM_STRING == "HA3")
      )
      {
        // return zero number of bonds
        return nr_bonds;
      }

      // determine the number of bonds if atom is an amide hydrogen
      if( ATOM_STRING == "H")
      {
        return 3;
      }

      // use the atom_pos char to determine how many bonds away the hydrogen atom is from the CB
      switch( atom_pos)
      {
        case 'A':
          nr_bonds = 2;
          break;
        case 'B':
          nr_bonds = 1;
          break;
        case 'G':
          nr_bonds = 2;
          break;
        case 'D':
          nr_bonds = 3;
          break;
        case 'E':
          nr_bonds = 4;
          break;
        case 'Z':
          nr_bonds = 5;
          break;
        case 'H':
          nr_bonds = 6;
          break;

        // side chain atom position not determined so return undefined
        default:
          return util::GetUndefined< size_t>();
      }

      // end
      return nr_bonds;
    }

    //! @brief returns whether a hydrogen atom type is in the backbone
    //! @return whether a hydrogen atom type is in the backbone
    bool StatisticNoe::IsBackboneProton( const biol::AtomType &ATOM_TYPE) const
    {
      return IsBackboneH( ATOM_TYPE) || IsBackboneHA( ATOM_TYPE);
    }

    //! @brief returns whether a hydrogen atom type is H
    //! @return whether a hydrogen atom type is H
    bool StatisticNoe::IsBackboneH( const biol::AtomType &ATOM_TYPE) const
    {
      return ATOM_TYPE == biol::GetAtomTypes().H;
    }

    //! @brief returns whether a hydrogen atom type is HA
    //! @return whether a hydrogen atom type is HA
    bool StatisticNoe::IsBackboneHA( const biol::AtomType &ATOM_TYPE) const
    {
      return ATOM_TYPE == biol::GetAtomTypes().HA ||
        ATOM_TYPE == biol::GetAtomTypes().HA2 ||
        ATOM_TYPE == biol::GetAtomTypes().HA3;
    }

  } // namespace app
} // namespace bcl
