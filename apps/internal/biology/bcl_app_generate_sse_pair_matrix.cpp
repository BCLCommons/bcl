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
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "assemble/bcl_assemble_sse_pair_template.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_triplet.h"

namespace bcl
{
  namespace app
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateSSEPairMatrix
    //! @brief This app is for creating pairwise RMSD tables for SSE Pair Orientations or Fold Templates.
    //!
    //! @author weinerbe, alexanns
    //! @date 11/14/09
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateSSEPairMatrix :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for specifying input file with distances between objects
      util::ShPtr< command::FlagInterface> m_InputFileFlag;

      //! flag for setting the binning information
      util::ShPtr< command::FlagStatic> m_BinsFlag;
      util::ShPtr< command::ParameterInterface> m_BinStartParameter;
      util::ShPtr< command::ParameterInterface> m_BinSizeParameter;
      util::ShPtr< command::ParameterInterface> m_BinNumberParameter;

      //! flag for creating final file with transformation matrices and weights
      util::ShPtr< command::FlagStatic> m_GenerateTransformationMatricesFlag;
      util::ShPtr< command::ParameterInterface> m_NodeCentersFilenameParameter;
      util::ShPtr< command::ParameterInterface> m_OutputFilenameParameter;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GenerateSSEPairMatrix();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      virtual GenerateSSEPairMatrix *Clone() const
      {
        return new GenerateSSEPairMatrix( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
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
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief GenerateTransformationMatrices creates final transformation matrix file
      //! @param SSE_PAIR_INFORMATION sse information read in from input file
      void GenerateTransformationMatrices
      (
        const storage::Map< size_t, util::ShPtr< assemble::SSEPairTemplate> > &SSE_PAIR_INFORMATION
      ) const;

      //! @brief generates tables of pairwise RMSDs for each bin
      //! @param SORTED_PAIRS map of bin and ShPtrVector of SSE pair templates to be used to calculate pairwise RMSDs
      //! @param PREFIX prefix to be placed in front of output files
      void GeneratePairwiseRMSDs
      (
        const storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > > &SORTED_PAIRS,
        const std::string &PREFIX
      ) const;

      //! @brief OutputRMSDTables takes a map of tables and outputs the tables to the filename specified
      //! @param TABLE_MAP the map which has all the rmsd information in it for each loop length bin
      //! @param OUTPUT_FILENAME the base filename which should be used to output each table in TABLE_MAP
      static void OutputRMSDTable
      (
        const storage::Table< double> &TABLE_MAP, const std::string &OUTPUT_FILENAME
      );

      //! @brief gets the coordinates from a SSE pair template
      //! @param PAIR_TEMPLATE SSE pair template containing the atom coordinates
      //! @return vector of coordinates from a SSE pair template
      static storage::Vector< linal::Vector3D> GetCoordinates( const assemble::SSEPairTemplate &PAIR_TEMPLATE);

      //! @brief gets a string of poly-As of a length determined by the SSType
      //! @param SS_TYPE SSType of SSE that the string will be used to create
      //! @return string of poly-As of a length determined by the SSType
      static std::string GetDummyString( const biol::SSType &SS_TYPE);

    private:

      static const ApplicationType GenerateSSEPairMatrix_Instance;

    }; // GenerateSSEPairMatrix

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> GenerateSSEPairMatrix::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input file flag
      sp_cmd->AddFlag( m_InputFileFlag);

      // binning information
      sp_cmd->AddFlag( m_BinsFlag);

      // flag for creating final file with transformation matrices and weights
      sp_cmd->AddFlag( m_GenerateTransformationMatricesFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateSSEPairMatrix::Main() const
    {
      // create input file stream "read"
      io::IFStream read;

      // create const string "input_filename" and initialize with the filename passed over the command line
      const std::string input_filename( m_InputFileFlag->GetFirstParameter()->GetValue());

      // open "input_filename" and bind it to "read"
      io::File::MustOpenIFStream( read, input_filename);

      // create the vector to hold the sse pair information
      storage::Vector< storage::Triplet< assemble::SSEGeometry, assemble::SSEGeometry, size_t> > sse_pairs;

      // fill "sse_pair_information" with the information from "input_filename"
      read >> sse_pairs;
      io::File::CloseClearFStream( read);

      // create maps for each of the contact types
      //           loop length                                        pair template      index
      storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > > sorted_helix_helix_pairs;
      storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > > sorted_helix_strand_pairs;
      storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > > sorted_strand_helix_pairs;
      storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > > sorted_strand_strand_pairs;

      // create map for all pairs, where the size_t is the index
      storage::Map< size_t, util::ShPtr< assemble::SSEPairTemplate> > all_pairs;

      // initialize the index
      size_t index( 0);

      // iterate through the sse pairs
      for
      (
        storage::Vector< storage::Triplet< assemble::SSEGeometry, assemble::SSEGeometry, size_t> >::const_iterator
          pair_itr( sse_pairs.Begin()), pair_itr_end( sse_pairs.End());
        pair_itr != pair_itr_end; ++pair_itr
      )
      {
        static const math::TransformationMatrix3D s_default_matrix;
        // if either sse geometry is undefined or the loop length is zero
        if
        (
          !pair_itr->First().IsDefined() || pair_itr->First().GetOrientation() == s_default_matrix ||
          !pair_itr->Second().IsDefined() || pair_itr->Second().GetOrientation() == s_default_matrix ||
          pair_itr->Third() == 0
        )
        {
          // skip to the next entry
          continue;
        }

        // create an SSEPairTemplate
        assemble::SSEPairTemplate pair_template
        (
          util::ShPtr< assemble::SSEGeometryInterface>( pair_itr->First().Clone()),
          util::ShPtr< assemble::SSEGeometryInterface>( pair_itr->Second().Clone()),
          pair_itr->Third()
        );

        // determine the current bin
        const size_t current_bin
        (
          ( double( pair_template.GetLoopLength()) - m_BinStartParameter->GetNumericalValue< double>()) /
          m_BinSizeParameter->GetNumericalValue< double>()
        );

        // if the current bin is not in the correct range
        if( current_bin > m_BinNumberParameter->GetNumericalValue< size_t>() - 1)
        {
          // skip to the next entry
          continue;
        }

        // construct the contact type
        contact::Type contact_type
        (
          contact::GetTypes().TypeFromSSTypes( *pair_template.GetFirstGeometry(), *pair_template.GetSecondGeometry())
        );

        // pushback the sse pair into the appropriate map
        if( contact_type == contact::GetTypes().HELIX_HELIX)
        {
          sorted_helix_helix_pairs[ current_bin].PushBack
          (
            util::ShPtr< storage::Pair< assemble::SSEPairTemplate, size_t> >
            (
              new storage::Pair< assemble::SSEPairTemplate, size_t>( pair_template, index)
            )
          );
        }
        else if( contact_type == contact::GetTypes().HELIX_SHEET)
        {
          sorted_helix_strand_pairs[ current_bin].PushBack
          (
            util::ShPtr< storage::Pair< assemble::SSEPairTemplate, size_t> >
            (
              new storage::Pair< assemble::SSEPairTemplate, size_t>( pair_template, index)
            )
          );
        }
        else if( contact_type == contact::GetTypes().SHEET_HELIX)
        {
          sorted_strand_helix_pairs[ current_bin].PushBack
          (
            util::ShPtr< storage::Pair< assemble::SSEPairTemplate, size_t> >
            (
              new storage::Pair< assemble::SSEPairTemplate, size_t>( pair_template, index)
            )
          );
        }
        else if( contact_type == contact::GetTypes().STRAND_STRAND)
        {
          sorted_strand_strand_pairs[ current_bin].PushBack
          (
            util::ShPtr< storage::Pair< assemble::SSEPairTemplate, size_t> >
            (
              new storage::Pair< assemble::SSEPairTemplate, size_t>( pair_template, index)
            )
          );
        }

        // if the generate transformation matrix flag was given
        if( m_GenerateTransformationMatricesFlag->GetFlag())
        {
          // also add to the total pairs
          all_pairs[ index] = util::ShPtr< assemble::SSEPairTemplate>( pair_template.Clone());
        }

        // increment the index
        ++index;
      }

      // if the generate transformation matrix flag was given
      if( m_GenerateTransformationMatricesFlag->GetFlag())
      {
        GenerateTransformationMatrices( all_pairs);
        return 0;
      }

      //  calculate the pairwise RMSDs and output the tables
      GeneratePairwiseRMSDs( sorted_helix_helix_pairs, "hh");
      GeneratePairwiseRMSDs( sorted_helix_strand_pairs, "hs");
      GeneratePairwiseRMSDs( sorted_strand_helix_pairs, "sh");
      GeneratePairwiseRMSDs( sorted_strand_strand_pairs, "ss");

      //successful end
      return 0;
    }

    //! @brief GenerateTransformationMatrices creates final transformation matrix file
    //! @param SSE_PAIR_INFORMATION sse information read in from input file
    void GenerateSSEPairMatrix::GenerateTransformationMatrices
    (
      const storage::Map< size_t, util::ShPtr< assemble::SSEPairTemplate> > &SSE_PAIR_INFORMATION
    ) const
    {
      // open input and output files
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_NodeCentersFilenameParameter->GetValue());
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_OutputFilenameParameter->GetValue());

      // iterate through the input file
      while( !read.eof())
      {
        // initialize temp variables
        size_t node_center;
        size_t weight;

        // read in node center and weight
        read >> node_center;
        read >> weight;

        // find the pair template in the map that matches the node center
        storage::Map< size_t, util::ShPtr< assemble::SSEPairTemplate> >::const_iterator find_itr
        (
          SSE_PAIR_INFORMATION.Find( node_center)
        );

        // the node center must be in the map, so assert that this is true
        BCL_Assert
        (
          find_itr != SSE_PAIR_INFORMATION.End(),
          "Node center of index, " + util::Format()( node_center) + " was not found in the list of acceptable indices"
        );

        // write out the pair information
        write << find_itr->second->GetPacking().GetContactType() << '\n';
        write << find_itr->second->GetLoopLength() << '\n';
        write << weight << '\n';

        // get a transformation matrix that holds the position of the second geometry whe the first is at the origin
        math::TransformationMatrix3D geometry_b( find_itr->second->GetSecondGeometry()->GetOrientation());
        geometry_b( math::Inverse( find_itr->second->GetFirstGeometry()->GetOrientation()));
        write << geometry_b << '\n';
      }
      io::File::CloseClearFStream( read);
      io::File::CloseClearFStream( write);
    }

    //! @brief generates tables of pairwise RMSDs for each bin
    //! @param SORTED_PAIRS map of bin and ShPtrVector of SSE pair templates to be used to calculate pairwise RMSDs
    //! @param PREFIX prefix to be placed in front of output files
    void GenerateSSEPairMatrix::GeneratePairwiseRMSDs
    (
      const storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > > &SORTED_PAIRS,
      const std::string &PREFIX
    ) const
    {
      // iterate through the map
      for
      (
        storage::Map< size_t, util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> > >::const_iterator
          map_itr( SORTED_PAIRS.Begin()), map_itr_end( SORTED_PAIRS.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // create vector of vector doubles "rows_and_values" which will hold all of the pairwise rmsd values calculated
        // for what will be the currentrmsd table
        storage::Vector< storage::Vector< double> > rows_and_values;

        // create vector "row_names" which will hold the string names of all of the sse pair fragments involved in the
        // current set of RMSD calculations
        storage::Vector< std::string> row_names;

        // iterate through the vector
        for
        (
          util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> >::const_iterator
            pair_itr( map_itr->second.Begin()), pair_itr_end( map_itr->second.End());
          pair_itr != pair_itr_end; ++pair_itr
        )
        {
          // add a vector to "rows_and_values" to hold the rmsds for the current row
          rows_and_values.PushBack( storage::Vector< double>());

          // add row name that is based on the line number from the input file
          row_names.PushBack( util::Format()( ( *pair_itr)->Second()));

          // store the coordinates from both sses
          storage::Vector< linal::Vector3D> pair_a_coords( GetCoordinates( ( *pair_itr)->First()));

          // iterate through the sse pairs further in order to get another pair for RMSD calculations
          for
          (
            util::ShPtrVector< storage::Pair< assemble::SSEPairTemplate, size_t> >::const_iterator
              pair_itr_b( map_itr->second.Begin());
            pair_itr_b != pair_itr_end; ++pair_itr_b
          )
          {
            // create variable to store rmsd value
            double calculated_rmsd;

            // create a table w/ and upper triangle, statement is true if the iterators are in the upper half
            if( pair_itr_b > pair_itr)
            {
              // store the coordinates
              storage::Vector< linal::Vector3D> pair_b_coords( GetCoordinates( ( *pair_itr_b)->First()));

              // create "rmsd_calculator" to do rmsd calculations
              const quality::RMSD rmsd_calculator;

              // calculate the rmsd
              calculated_rmsd = rmsd_calculator.CalculateMeasure
                (
                  util::ConvertToConstSiPtrVector< linal::Vector3D>( pair_a_coords),
                  util::ConvertToConstSiPtrVector< linal::Vector3D>( pair_b_coords)
                );
            }
            else
            {
              // iterators are in lower triangle, so store 0.0
              calculated_rmsd = 0.0;
            }

            // calculate and pushback to the vector of RMSDs that is being created for the current row
            rows_and_values.ReverseBegin()->PushBack( calculated_rmsd);
          }
        }

        // create TableHeader "table_header" and initialize with the vector of row names "row_names" which was compiled
        // all of the names of sse fragment pairs that were involved in the set of RMSD calculations
        storage::TableHeader table_header( row_names);

        // create Table "current_rmsd_table" and initialize with "table_header"
        storage::Table< double> current_rmsd_table( table_header);

        // create const_iterator "row_names_itr" and "row_names_itr_end" to iterate through the row names and help
        // to fill the rows of "current_rmsd_table" in the for loop immediately below
        storage::Vector< std::string>::const_iterator row_names_itr( row_names.Begin()), row_names_itr_end( row_names.End());

        // fill "current_rmsd_table" with the pairwise rmsd calculations that are contained in "rows_and_values"
        for
        (
          storage::Vector< storage::Vector< double> >::const_iterator
            row_itr( rows_and_values.Begin()), row_itr_end( rows_and_values.End());
          row_itr != row_itr_end && row_names_itr != row_names_itr_end;
          ++row_itr, ++row_names_itr
        )
        {
          // insert the current row of rmsd values into "current_rmsd_table"
          current_rmsd_table.InsertRow( *row_names_itr, *row_itr, true);
        }

        // append bin number to the given filename
        const std::string current_filename( PREFIX + "_" + util::Format()( map_itr->first) + ".out");

        // output the table using OUTPUT_FILENAME
        OutputRMSDTable( current_rmsd_table, current_filename);
      }
    }

    //! @brief OutputRMSDTables takes a map of tables and outputs the tables to the filename specified
    //! @param TABLE_MAP the map which has all the rmsd information in it for each loop length bin
    //! @param OUTPUT_FILENAME the base filename which should be used to output each table in TABLE_MAP
    void GenerateSSEPairMatrix::OutputRMSDTable
    (
      const storage::Table< double> &TABLE_MAP, const std::string &OUTPUT_FILENAME
    )
    {
      // open output file stream and write the table
      io::OFStream write;
      io::File::MustOpenOFStream( write, OUTPUT_FILENAME);
      TABLE_MAP.WriteFormatted( write);
      io::File::CloseClearFStream( write);
    }

    //! @brief gets the coordinates from a SSE pair template
    //! @param PAIR_TEMPLATE SSE pair template containing the atom coordinates
    //! @return vector of coordinates from a SSE pair template
    storage::Vector< linal::Vector3D> GenerateSSEPairMatrix::GetCoordinates
    (
      const assemble::SSEPairTemplate &PAIR_TEMPLATE
    )
    {
      // create dummy sses for the geometries
      biol::AASequence dummy_sequence_a
      (
        biol::AASequenceFactory::BuildSequenceFromFASTAString
        (
          GetDummyString( PAIR_TEMPLATE.GetFirstGeometry()->GetType()),
          biol::GetAAClasses().e_AABackBone
        )
      );
      assemble::SSE dummy_sse_a( dummy_sequence_a, PAIR_TEMPLATE.GetFirstGeometry()->GetType());
      biol::AASequence dummy_sequence_b
      (
        biol::AASequenceFactory::BuildSequenceFromFASTAString
        (
          GetDummyString( PAIR_TEMPLATE.GetSecondGeometry()->GetType()),
          biol::GetAAClasses().e_AABackBone
        )
      );
      assemble::SSE dummy_sse_b( dummy_sequence_b, PAIR_TEMPLATE.GetSecondGeometry()->GetType());

      // idealize the sses then move them to the positions given by the geometries
      dummy_sse_a.SetToIdealConformationAtOrigin();
      dummy_sse_b.SetToIdealConformationAtOrigin();
      dummy_sse_a.Transform( PAIR_TEMPLATE.GetFirstGeometry()->GetOrientation());
      dummy_sse_b.Transform( PAIR_TEMPLATE.GetSecondGeometry()->GetOrientation());

      // store the coordinates from both sses
      util::SiPtrVector< const linal::Vector3D> sp_pair_coords( dummy_sse_a.GetAtomCoordinates());
      sp_pair_coords.Append( dummy_sse_b.GetAtomCoordinates());
      storage::Vector< linal::Vector3D> pair_coords;

      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator coord_itr( sp_pair_coords.Begin()),
          coord_itr_end( sp_pair_coords.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        pair_coords.PushBack( **coord_itr);
      }

      // end
      return pair_coords;
    }

    //! @brief gets a string of poly-As of a length determined by the SSType
    //! @param SS_TYPE SSType of SSE that the string will be used to create
    //! @return string of poly-As of a length determined by the SSType
    std::string GenerateSSEPairMatrix::GetDummyString( const biol::SSType &SS_TYPE)
    {
      if( SS_TYPE == biol::GetSSTypes().HELIX)
      {
        return "AAAAA";
      }
      else if( SS_TYPE == biol::GetSSTypes().STRAND)
      {
        return "AAA";
      }

      return "";
    }

    // default constructor
    GenerateSSEPairMatrix::GenerateSSEPairMatrix() :
      m_InputFileFlag
      (
        new command::FlagStatic
        (
          "input_file", "sse pair or fold template information list",
          command::Parameter
          (
            "input_filename",
            "filename of sse pair or fold template information list"
          )
        )
      ),
      m_BinsFlag
      (
        new command::FlagStatic
        (
          "bins",
          "flag for specifying the binning to be used"
        )
      ),
      m_BinStartParameter
      (
        new command::Parameter
        (
          "loop_start_bin",
          "This is the starting loop length that will be used",
          "1"
        )
      ),
      m_BinSizeParameter
      (
        new command::Parameter
        (
          "loop_bin_size",
          "This is the bin size that will be used",
          "1"
        )
      ),
      m_BinNumberParameter
      (
        new command::Parameter
        (
          "loop_number_bins",
          "This is the number of bins that the loop length will be split into",
          "20"
        )
      ),
      m_GenerateTransformationMatricesFlag
      (
        new command::FlagStatic
        (
          "generate_transformation_matrices", "generate file with tranformation matrices and weights"
        )
      ),
      m_NodeCentersFilenameParameter
      (
        new command::Parameter( "node_centers_filename", "file containing node centers to be used", "nodes.list")
      ),
      m_OutputFilenameParameter
      (
        new command::Parameter
        (
          "output_filename", "output file containing the transformation matrices", "sse_pair_matrices.out"
        )
      )
    {
      // attach parameters to flags
      m_BinsFlag->PushBack( m_BinStartParameter);
      m_BinsFlag->PushBack( m_BinSizeParameter);
      m_BinsFlag->PushBack( m_BinNumberParameter);
      m_GenerateTransformationMatricesFlag->PushBack( m_NodeCentersFilenameParameter);
      m_GenerateTransformationMatricesFlag->PushBack( m_OutputFilenameParameter);
    }

    const ApplicationType GenerateSSEPairMatrix::GenerateSSEPairMatrix_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateSSEPairMatrix(), GetAppGroups().e_Protein)
    );

  } // namespace app
} // namespace bcl

