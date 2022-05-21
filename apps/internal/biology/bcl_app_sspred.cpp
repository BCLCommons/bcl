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
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"

namespace bcl
{
  namespace app
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSPred
    //! @brief Calculates statistic for ss prediction methods
    //! @details Read, analyze and generate statistics for any secondary structure prediction for the given sequences
    //!
    //! @author karakam
    //! @date 09/08/08
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSPred :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! input for a single file
      util::ShPtr< command::FlagStatic> m_FastaFlag;
      util::ShPtr< command::ParameterInterface> m_FastaParam;
      util::ShPtr< command::ParameterInterface> m_FastaOligomericStateParam;

      //! input for multiple files
      util::ShPtr< command::FlagStatic> m_FastaListFlag;
      util::ShPtr< command::ParameterInterface> m_FastaListParam;

      //! path for input files
      util::ShPtr< command::FlagStatic> m_PathFlag;
      util::ShPtr< command::ParameterInterface> m_PathParam;

      //! flag for using pdb path hierarchy for all input files
      util::ShPtr< command::FlagInterface> m_PdbPathHierarchy;

      //! path for output files
      util::ShPtr< command::FlagStatic> m_OutputPathFlag;
      util::ShPtr< command::ParameterInterface> m_OutputPathParam;

      //! flag for calculating statistics
      util::ShPtr< command::FlagInterface> m_StatisticsFlag;

      //! flag for reading predictions from different SSMethods
      util::ShPtr< command::FlagStatic> m_PdbListFlag;
      util::ShPtr< command::ParameterInterface> m_PdbListParam;

      //! flag for providing a oligomeric state dictionary
      util::ShPtr< command::FlagStatic> m_OligomericStateDictionaryFlag;
      util::ShPtr< command::ParameterInterface> m_OligomericStateDictionaryParam;

      //! flag to set the path where membrane related data can be found
      util::ShPtr< command::FlagStatic> m_MembraneDataPathFlag;
      util::ShPtr< command::ParameterInterface> m_MembraneDataPathParam;

      //! flag to use pdbtmxml files for transformation matrices
      util::ShPtr< command::FlagInterface> m_PdbtmXmlFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      SSPred();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      SSPred *Clone() const
      {
        return new SSPred( *this);
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

      //! @brief read all pdb tags and three membrane thickness
      //! @return list of triplets containing a string for each pdb filename, and a membrane object and a transformation
      storage::List< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
      ReadPdbTagsAndMembranes() const;

      //! @brief initializes the command object for that executable
      //! @return initalized command object
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief Main function
      //! @return return value of the application
      int Main() const;

      //! @brief calculates percentages by normalizing over rows
      //! @param TABLE matrix to be used
      //! @return new matrix with percentages
      linal::Matrix< double> CalculatePercentagesByRows( const linal::Matrix< double> &TABLE) const;

      //! @brief calculates percentages by normalizing over columns
      //! @param TABLE matrix to be used
      //! @return new matrix with percentages
      linal::Matrix< double> CalculatePercentagesByColumns( const linal::Matrix< double> &TABLE) const;

      //! @brief outputs a 9 state table
      //! @param NINE_STATE_TABLE 9 state table
      //! @param FORMAT formatting to be used
      void Print9StateTable( const linal::Matrix< double> &NINE_STATE_TABLE, const util::Format &FORMAT) const;

      //! @brief outputs a 3 state sspred table
      //! @param THREE_STATE_TABLE 3 state table
      //! @param FORMAT formatting to be used
      void Print3StateSSTable( const linal::Matrix< double> &THREE_STATE_TABLE, const util::Format &FORMAT) const;

      //! @brief outputs a 3 state tm table
      //! @param THREE_STATE_TABLE 3 state table
      //! @param FORMAT formatting to be used
      void Print3StateTMTable( const linal::Matrix< double> &THREE_STATE_TABLE, const util::Format &FORMAT) const;

      //! @brief outputs a 2 state table
      //! @param TWO_STATE_TABLE 2 state table
      //! @param STATE_NAME name of the state of interest
      //! @param FORMAT formatting to be used
      void Print2StateTable( const linal::Matrix< double> &TWO_STATE_TABLE, const std::string &STATE_NAME, const util::Format &FORMAT) const;

      //! @brief calculates percentage of correct predictions
      //! @param TABLE Matrix to be processed
      //! @return percentage of correct predictions
      double CalculatePercentageCorrectPredictions( const linal::Matrix< double> &TABLE) const;

      //! @brief calculates information gain for a table
      //! @param TABLE Matrix to be processed
      //! @return information gain
      double CalculateInformationGain( const linal::Matrix< double> &TABLE) const;

      //! @brief calculates mutual information for a table
      //! @param TABLE Matrix to be processed
      //! @return mutual information
      double CalculateMutualInformation( const linal::Matrix< double> &TABLE) const;

      //! @brief collapses table to a 2 state table
      //! @param TABLE table of interest
      //! @param INDEX index of the state of interest
      //! @return collapsed 2state table
      linal::Matrix< double> CollapseTableTo2State( const linal::Matrix< double> &TABLE, const size_t INDEX) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      //! static instance of this class
      static const ApplicationType SSPred_Instance;

    }; // class SSPred

    //! @brief Main function
    //! @return return value of the application
    int SSPred::Main() const
    {
      // initialize read and write streams
      io::IFStream read;
      io::OFStream write;

      // initialize pdb factory and handler
      pdb::Factory pdb_factory( biol::GetAAClasses().e_AABackBone);

      // initialize sse_min_sizes
      storage::Map< biol::SSType, size_t> sse_min_sizes;
      sse_min_sizes[ biol::GetSSTypes().HELIX] = 0;
      sse_min_sizes[ biol::GetSSTypes().STRAND] = 0;
      sse_min_sizes[ biol::GetSSTypes().COIL] = 0;

      // initialize the set of SSMethods to be used
      const storage::Set< sspred::Method> ss_methods( m_StatisticsFlag->GetObjectSet< sspred::Method>());

      // initialize list to store pdb names and membrane associated information
      storage::List< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
      pdb_name_membrane_list( ReadPdbTagsAndMembranes());

      // initialize matrices to store accuracy averages and prediction counts
      storage::Vector< linal::Matrix< double> > total_prediction_counts
      (
        sspred::GetMethods().GetEnumCount(), linal::Matrix< double>( 9, 9, double( 0.0))
      );
      storage::Vector< linal::Matrix< double> > total_prediction_counts_sspred
      (
        sspred::GetMethods().GetEnumCount(), linal::Matrix< double>( 3, 3, double( 0.0))
      );
      storage::Vector< linal::Matrix< double> > total_prediction_counts_tmpred
      (
        sspred::GetMethods().GetEnumCount(), linal::Matrix< double>( 3, 3, double( 0.0))
      );

      // iterate over the pdb list
      for
      (
        storage::List< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >::const_iterator
          pdb_itr( pdb_name_membrane_list.Begin()), pdb_itr_end( pdb_name_membrane_list.End());
        pdb_itr != pdb_itr_end; ++pdb_itr
      )
      {
        util::GetLogger() << "PDB: " << pdb_itr->First() << '\n';

        // initialize path
        std::string input_path( m_PathParam->GetValue() + PATH_SEPARATOR);
        // if pdb hiearchy is used
        if( m_PdbPathHierarchy->GetFlag())
        {
          input_path += pdb_itr->First().substr( 1, 2) + PATH_SEPARATOR;
        }

        // read pdb
        io::File::MustOpenIFStream( read, input_path + pdb_itr->First() + ".pdb");
        std::string prefix_without_chainid( pdb_itr->First().substr( 0, 4));
        pdb::Handler pdb( read);
        assemble::ProteinModel this_model( pdb_factory.ProteinModelFromPDB( pdb, sse_min_sizes));
        io::File::CloseClearFStream( read);

        // read all SSMethods
        sspred::MethodHandler::ReadPredictionsForProteinModel
        (
          ss_methods, this_model, prefix_without_chainid, input_path
        );

        // transform model if third column in input file is transformation matrix
        this_model.Transform( pdb_itr->Third());

        // set membrane
        util::ShPtr< assemble::ProteinModelData> sp_data( this_model.GetProteinModelData());
        sp_data->Insert( assemble::ProteinModelData::e_Membrane, pdb_itr->Second());
        this_model.SetProteinModelData( sp_data);
        BCL_MessageStd( util::Format()( pdb_itr->Second()));

        // update PDB "prediction" w/ environment types
        sspred::PDB::SetEnvironmentTypes( this_model);

        // iterate over methods
        for
        (
          storage::Set< sspred::Method>::const_iterator method_itr( ss_methods.Begin()), method_itr_end( ss_methods.End());
          method_itr != method_itr_end; ++method_itr
        )
        {
          util::GetLogger() << "METHOD: " << method_itr->GetName() << '\n';

          // initialize counts for this model
          linal::Matrix< double> model_prediction_counts( 9, 9, double( 0.0));
          linal::Matrix< double> model_prediction_counts_sspred( 3, 3, double( 0.0));
          linal::Matrix< double> model_prediction_counts_tmpred( 3, 3, double( 0.0));

          // get all SSEs
          const util::SiPtrVector< const assemble::SSE> all_sses( this_model.GetSSEs());
          io::File::CloseClearFStream( read);

          // iterate over all SSEs
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
            sse_itr != sse_itr_end; ++sse_itr
          )
          {
            // store sstype
            const biol::SSType this_sstype( ( *sse_itr)->GetType());

            // iterate over all residues
            for
            (
              biol::AASequence::const_iterator aa_itr( ( *sse_itr)->GetData().Begin()),
                aa_itr_end( ( *sse_itr)->GetData().End());
              aa_itr != aa_itr_end; ++aa_itr
            )
            {
              // if CA coordinates are not defined
              if( !( *aa_itr)->GetCA().GetCoordinates().IsDefined())
              {
                // then skip this amino acid
                continue;
              }

              // determine environment type
              biol::EnvironmentType environment
              (
                pdb_itr->Second()->DetermineEnvironmentType( ( *aa_itr)->GetCA().GetCoordinates())
              );

              // determine one state prediction
              const storage::Pair< biol::SSType, biol::EnvironmentType> prediction
              (
                ( *aa_itr)->GetSSPrediction( *method_itr)->GetOneStateSSTMPrediction()
              );

              BCL_MessageDbg
              (
                "AA: " + ( *aa_itr)->GetIdentification() + " " +
                "ORIG: " + environment.GetName() + " " + this_sstype.GetName() + " " +
                "PRED: " + prediction.Second().GetName() + " " + prediction.First().GetName()
              );

              // sum up the counts
              model_prediction_counts( ( environment->GetReducedIndex() * 3 + this_sstype), ( prediction.Second()->GetReducedIndex() * 3 + prediction.First())) += 1;
              model_prediction_counts_sspred( this_sstype, prediction.First()) += 1;
              model_prediction_counts_tmpred( environment->GetReducedIndex(), prediction.Second()->GetReducedIndex()) += 1;
            }
          }

          // print the tables and percentage tables
          Print9StateTable( model_prediction_counts, util::Format().W( 4));
          util::GetLogger() << "\n=================================================================\n\n";
          Print9StateTable( CalculatePercentagesByRows( model_prediction_counts), util::Format().W( 6).FFP( 2));
          util::GetLogger() << "\n=================================================================\n\n";
          Print3StateSSTable( model_prediction_counts_sspred, util::Format().W( 4));
          util::GetLogger() << "\n=================================================================\n\n";
          Print3StateSSTable( CalculatePercentagesByRows( model_prediction_counts_sspred), util::Format().W( 6).FFP( 2));
          util::GetLogger() << "\n=================================================================\n\n";
          Print3StateTMTable( model_prediction_counts_tmpred, util::Format().W( 4));
          util::GetLogger() << "\n=================================================================\n\n";
          Print3StateTMTable( CalculatePercentagesByRows( model_prediction_counts_tmpred), util::Format().W( 6).FFP( 2));
          util::GetLogger() << "\n=================================================================\n\n";
          // collapse to transmembrane helix
          Print2StateTable
          (
            CollapseTableTo2State( model_prediction_counts, biol::GetSSTypes().HELIX + 3 * biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex()),
            "MC-H",
            util::Format().W( 4)
          );
          util::GetLogger() << "\n=================================================================\n\n";
          // collapse to transmembrane strand
          Print2StateTable
          (
            CollapseTableTo2State( model_prediction_counts, biol::GetSSTypes().STRAND + 3 * biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex()),
            "MC-E",
            util::Format().W( 4)
          );
          util::GetLogger() << "\n=================================================================\n\n";

          // update the counts in the global matrices
          total_prediction_counts( *method_itr) += model_prediction_counts;
          total_prediction_counts_sspred( *method_itr) += model_prediction_counts_sspred;
          total_prediction_counts_tmpred( *method_itr) += model_prediction_counts_tmpred;

        } // all methods
      } // pdb_itr

      util::GetLogger() << "THE TOTAL COUNTS AND PERCENTAGES" << '\n';

      // iterate over methods
      for
      (
        storage::Set< sspred::Method>::const_iterator method_itr( ss_methods.Begin()), method_itr_end( ss_methods.End());
        method_itr != method_itr_end; ++method_itr
      )
      {
        util::GetLogger() << "METHOD: " << method_itr->GetName() << '\n';
        // print the tables
        util::GetLogger() << "\n=================================================================\n\n";
        Print9StateTable( total_prediction_counts( *method_itr), util::Format().W( 7));
        util::GetLogger() << "\n=================================================================\n\n";
        Print9StateTable( CalculatePercentagesByRows( total_prediction_counts( *method_itr)), util::Format().W( 6).FFP( 2));
        util::GetLogger() << "\n=================================================================\n\n";
        Print3StateSSTable( total_prediction_counts_sspred( *method_itr), util::Format().W( 7));
        util::GetLogger() << "\n=================================================================\n\n";
        Print3StateSSTable( CalculatePercentagesByRows( total_prediction_counts_sspred( *method_itr)), util::Format().W( 6).FFP( 2));
        util::GetLogger() << "\n=================================================================\n\n";
        Print3StateTMTable( total_prediction_counts_tmpred( *method_itr), util::Format().W( 7));
        util::GetLogger() << "\n=================================================================\n\n";
        Print3StateTMTable( CalculatePercentagesByRows( total_prediction_counts_tmpred( *method_itr)), util::Format().W( 6).FFP( 2));
        util::GetLogger() << "\n=================================================================\n\n";
        // collapse to transmembrane helix
        Print2StateTable
        (
          CollapseTableTo2State( total_prediction_counts( *method_itr), biol::GetSSTypes().HELIX + 3 * biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex()),
          "MC-H",
          util::Format().W( 7)
        );
        util::GetLogger() << "\n=================================================================\n\n";
        // collapse to transmembrane strand
        Print2StateTable
        (
          CollapseTableTo2State( total_prediction_counts( *method_itr), biol::GetSSTypes().STRAND + 3 * biol::GetEnvironmentTypes().e_MembraneCore->GetReducedIndex()),
          "MC-E",
          util::Format().W( 7)
        );
        util::GetLogger() << "\n=================================================================\n\n";
      }

      // iterate over methods
      for
      (
        storage::Set< sspred::Method>::const_iterator method_itr( ss_methods.Begin()), method_itr_end( ss_methods.End());
        method_itr != method_itr_end; ++method_itr
      )
      {
        util::GetLogger() << "METHOD: " << method_itr->GetName() << '\n';
        // print the tables
        util::GetLogger() << "\t9 state information gain: " << CalculateInformationGain( total_prediction_counts( *method_itr)) << '\n';
        util::GetLogger() << "\t9 state mutual information: " << CalculateMutualInformation( total_prediction_counts( *method_itr)) << '\n';
        util::GetLogger() << "\tSSPred  information gain: " << CalculateInformationGain( total_prediction_counts_sspred( *method_itr)) << '\n';
        util::GetLogger() << "\tSSPred  mutual information: " << CalculateMutualInformation( total_prediction_counts_sspred( *method_itr)) << '\n';
        util::GetLogger() << "\tTMPred  information gain: " << CalculateInformationGain( total_prediction_counts_tmpred( *method_itr)) << '\n';
        util::GetLogger() << "\tTMPred  mutual information: " << CalculateMutualInformation( total_prediction_counts_tmpred( *method_itr)) << '\n';
      }

      return 0;
    }

    //! @brief read all pdb tags and three membrane thickness
    //! @return list of triplets containing a string for each pdb filename, and a membrane object and a transformation
    storage::List< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
    SSPred::ReadPdbTagsAndMembranes() const
    {
      // store membrane normals in storage vector of a storage::Pair< std::string, double> of pdb-name and its double value
      // for soluble proteins the membrane just does not have a core nor a transition regions or gaps
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_PdbListParam->GetValue());

      // store the .pdb tag as a string
      // store the membrane
      // store the transformation
      storage::List< storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> >
        pdb_membrane_transformation_list;

      // read line by line
      std::string this_line;
      while( std::getline( read, this_line) && !read.eof())
      {
        // initialize triplet
        storage::Triplet< std::string, util::ShPtr< biol::Membrane>, math::TransformationMatrix3D> name_membrane_trans;

        // split line
        storage::Vector< std::string> string_list( util::SplitString( this_line));

        // get the pdb name
        name_membrane_trans.First() = string_list( 0);

        // sometimes read.eof() does not return false, even if end is reached
        // so it is necessary to check if a name was actually read
        if( name_membrane_trans.First().empty())
        {
          break;
        }

        // if string has only one member or soluble is specified then it is a soluble protein
        if( string_list.GetSize() == 1 || string_list( 1) == "Soluble")
        {
          name_membrane_trans.Second() = util::ShPtr< biol::Membrane>( new biol::Membrane( 0.0, 0.0, 0.0));
        }
        // if membrane proteins
        else if( string_list( 1) == "Membrane")
        {
          // if there are 3 strings and pdbtm xml flag is given
          if( string_list.GetSize() == 2 && m_PdbtmXmlFlag->GetFlag())
          {
            io::IFStream read_xml;
            io::File::MustOpenIFStream( read_xml, m_MembraneDataPathParam->GetValue() + PATH_SEPARATOR + string_list( 2));
            const storage::Pair< biol::Membrane, math::TransformationMatrix3D>
              membrane_transformation
              (
                biol::Membrane::MembraneAndTransformationFromPDBTMXML
                (
                  read_xml,
                  biol::Membrane::GetParameterTransitionThickness()->GetNumericalValue< double>(),
                  biol::Membrane::GetParameterGapThickness()->GetNumericalValue< double>()
                )
              );
            name_membrane_trans.Second() = util::ShPtr< biol::Membrane>
            (
              membrane_transformation.First().Clone()
            );
            name_membrane_trans.Third() = membrane_transformation.Second();
          }
          // else membrane protein but no pdbtm xml files is given is given
          else
          {
            name_membrane_trans.Second() = util::ShPtr< biol::Membrane>
            (
              new biol::Membrane( biol::Membrane::GetCommandLineMembrane())
            );
          }
        }
        else
        {
          BCL_Exit( "each line in pdb list should have 1 (for solubles/membrane) or 2 (for membrane with xml) entries", -1);
        }

        // store current name and membrane
        pdb_membrane_transformation_list.PushBack( name_membrane_trans);
      }
      io::File::CloseClearFStream( read);

      // end
      return pdb_membrane_transformation_list;
    }

    //! @brief calculates percentages by normalizing over rows
    //! @param TABLE matrix to be used
    //! @return new matrix with percentages
    linal::Matrix< double> SSPred::CalculatePercentagesByRows( const linal::Matrix< double> &TABLE) const
    {
      // initialize matrix to be returned
      linal::Matrix< double> this_matrix( TABLE);

      // iterate over the rows
      for( size_t row( 0); row < TABLE.GetNumberRows(); ++row)
      {
        // initialize sum
        double this_sum( 0);

        // iterate over rows to find the sum
        for( size_t column( 0); column < TABLE.GetNumberCols(); ++column)
        {
          this_sum += TABLE( row, column);
        }
        // if the sum is equal to 0 then skip
        if( this_sum == 0)
        {
          continue;
        }
        // otherwise iterate again to normalize
        for( size_t column( 0); column < TABLE.GetNumberCols(); ++column)
        {
          this_matrix( row, column) /= this_sum / 100;
        }
      }
      // end
      return this_matrix;
    }

    //! @brief calculates percentages by normalizing over columns
    //! @param TABLE matrix to be used
    //! @return new matrix with percentages
    linal::Matrix< double> SSPred::CalculatePercentagesByColumns( const linal::Matrix< double> &TABLE) const
    {
      // initialize matrix to be returned
      linal::Matrix< double> this_matrix( TABLE);

      // iterate over the rows
      for( size_t col( 0); col < TABLE.GetNumberCols(); ++col)
      {
        // initialize sum
        double this_sum( 0);

        // iterate over rows to find the sum
        for( size_t row( 0); row < TABLE.GetNumberRows(); ++row)
        {
          this_sum += TABLE( row, col);
        }
        // iterate again to normalize
        for( size_t row( 0); row < TABLE.GetNumberRows(); ++row)
        {
          this_matrix = TABLE( row, col) / this_sum;
        }
      }
      // end
      return this_matrix;
    }

    //! @brief outputs a 9 state table
    //! @param NINE_STATE_TABLE 9 state table
    void SSPred::Print9StateTable( const linal::Matrix< double> &NINE_STATE_TABLE, const util::Format &FORMAT) const
    {
      util::GetLogger() << "9state    ";
      for( size_t i( 0); i < 3; ++i)
      {
        for( size_t j( 0); j < 3; ++j)
        {
          util::GetLogger() << biol::GetEnvironmentTypes().GetReducedTypes()( i)->GetTwoLetterCode() << "-" << biol::SSType( j)->GetOneLetterCode() << '\t';
        }
      }
      util::GetLogger() << '\n';
      for( size_t i( 0); i < 3; ++i)
      {
        for( size_t j( 0); j < 3; ++j)
        {
          // calculate and store row number
          const size_t row_no( i * 3 + j);
          util::GetLogger() << biol::GetEnvironmentTypes().GetReducedTypes()( i)->GetTwoLetterCode() << "-" << biol::SSType( j)->GetOneLetterCode() << '\t';

          for( size_t k( 0); k < 3; ++k)
          {
            for( size_t l( 0); l < 3; ++l)
            {
              util::GetLogger() << FORMAT( NINE_STATE_TABLE( row_no, k * 3 + l)) << '\t';
            }
          }
          util::GetLogger() << '\n';
        }
      }
      util::GetLogger() << '\n';
    }

    //! @brief outputs a 3 state sspred table
    //! @param THREE_STATE_TABLE 3 state table
    void SSPred::Print3StateSSTable( const linal::Matrix< double> &THREE_STATE_TABLE, const util::Format &FORMAT) const
    {
      util::GetLogger() << "sspred    ";
      for( size_t i( 0); i < 3; ++i)
      {
        util::GetLogger() << biol::SSType( i)->GetOneLetterCode() << '\t';
      }
      util::GetLogger() << '\n';
      for( size_t i( 0); i < 3; ++i)
      {
        util::GetLogger() << biol::SSType( i)->GetOneLetterCode() << '\t';
        for( size_t j( 0); j < 3; ++j)
        {
          util::GetLogger() << FORMAT( THREE_STATE_TABLE( i, j)) << '\t';
        }
        util::GetLogger() << '\n';
      }
      util::GetLogger() << '\n';
    }

    //! @brief outputs a 3 state tm table
    //! @param THREE_STATE_TABLE 3 state table
    void SSPred::Print3StateTMTable( const linal::Matrix< double> &THREE_STATE_TABLE, const util::Format &FORMAT) const
    {
      util::GetLogger() << "tmpred    ";
      for( size_t i( 0); i < 3; ++i)
      {
        util::GetLogger() << biol::GetEnvironmentTypes().GetReducedTypes()( i)->GetTwoLetterCode() << '\t';
      }
      util::GetLogger() << '\n';
      for( size_t i( 0); i < 3; ++i)
      {
        util::GetLogger() << biol::GetEnvironmentTypes().GetReducedTypes()( i)->GetTwoLetterCode() << '\t';
        for( size_t j( 0); j < 3; ++j)
        {
          util::GetLogger() << FORMAT( THREE_STATE_TABLE( i, j)) << '\t';
        }
        util::GetLogger() << '\n';
      }
      util::GetLogger() << '\n';
    }

    //! @brief outputs a 2 state table
    //! @param TWO_STATE_TABLE 2 state table
    //! @param STATE_NAME name of the state of interest
    void SSPred::Print2StateTable( const linal::Matrix< double> &TWO_STATE_TABLE, const std::string &STATE_NAME, const util::Format &FORMAT) const
    {
      util::GetLogger() << "TMPred    " << STATE_NAME << '\t' << "else" << '\n';
      util::GetLogger() << STATE_NAME << '\t' << TWO_STATE_TABLE( 0, 0) << '\t' << TWO_STATE_TABLE( 0, 1) << '\n';
      util::GetLogger() << "else\t"  << TWO_STATE_TABLE( 1, 0) << '\t' << TWO_STATE_TABLE( 1, 1) << '\n';
      util::GetLogger() << '\n';
    }

    //! @brief calculates percentage of correct predictions
    //! @param TABLE Matrix to be processed
    //! @return percentage of correct predictions
    double SSPred::CalculatePercentageCorrectPredictions( const linal::Matrix< double> &TABLE) const
    {
      return TABLE.Trace() * double( 100.0) / TABLE.Sum();
    }

    //! @brief calculates information gain for a table
    //! @param TABLE Matrix to be processed
    //! @return information gain
    double SSPred::CalculateInformationGain( const linal::Matrix< double> &TABLE) const
    {
      // store the background probability
      const double background_probability( 1.0 / double( TABLE.GetNumberOfElements()));

      // calculate the percentage table ( over rows)
      linal::Vector< double> diagonal( CalculatePercentagesByRows( TABLE).GetDiagonal());
      // add pseudocount to diagonal
      diagonal += std::pow( 10.0, -15.0);

      // initialize information gain
      double information_gain( 0);

      // iterate over the diagonal elements
      for( double *ptr( diagonal.Begin()), *ptr_end( diagonal.End()); ptr != ptr_end; ++ptr)
      {
        information_gain += std::log( *ptr / ( double( 100.0) * background_probability));
      }

      // normalize by the size of diagonal
      information_gain /= diagonal.GetSize();
      information_gain /= std::log( 2.0);

      // end
      return information_gain;
    }

    //! @brief calculates mutual information for a table
    //! @param TABLE Matrix to be processed
    //! @return mutual information
    double SSPred::CalculateMutualInformation( const linal::Matrix< double> &TABLE) const
    {
      // make a copy of the table
      linal::Matrix< double> this_table( TABLE);
      this_table += std::pow( 10.0, -15.0);

      // initialize mutual information
      double mutual_information( 0);

      // calculate percentages by row and column
      const linal::Matrix< double> percentages_rows( CalculatePercentagesByRows( this_table));
      const linal::Matrix< double> percentages_cols( CalculatePercentagesByColumns( this_table));

      // store the sum
      const double sum( this_table.Sum());

      // iterate over the rows and columns
      for( size_t row( 0); row < this_table.GetNumberRows(); ++row)
      {
        for( size_t col( 0); col < this_table.GetNumberCols(); ++col)
        {
          // sum over the mutual information for each element in the matrix
          mutual_information +=
            std::log( ( this_table( row, col) / sum) / ( percentages_rows( row, col) * percentages_cols( row, col) / 10000.0));
        }
      }
      // normalize by size
      mutual_information /= this_table.GetNumberOfElements();
      mutual_information /= std::log( 2.0);

      // end
      return mutual_information;
    }

    //! @brief collapses table to a 2 state table
    //! @param TABLE table of interest
    //! @param INDEX index of the state of interest
    //! @return collapsed 2state table
    linal::Matrix< double> SSPred::CollapseTableTo2State( const linal::Matrix< double> &TABLE, const size_t INDEX) const
    {
      // initialize matrix to be returned
      linal::Matrix< double> two_state_table( 2, 2, double( 0.0));

      // store the count of correct state
      const size_t count_correct_state( TABLE( INDEX, INDEX));

      // set the values
      two_state_table( 0, 0) = count_correct_state;
      two_state_table( 1, 0) = TABLE.GetCol( INDEX).Sum() - count_correct_state;
      two_state_table( 0, 1) = TABLE.GetRow( INDEX).Sum() - count_correct_state;
      // since 1,1 is set to 0 in the initializer two_state_table.Sum() gives the total count for the remainder of the TABLE
      two_state_table( 1, 1) = TABLE.Sum() - two_state_table.Sum();

      // end
      return two_state_table;
    }

    //! @brief initializes the command object for that executable
    //! @return initialized command object
    util::ShPtr< command::Command> SSPred::InitializeCommand() const
    {
      // initialize a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());
      // input for a single file
      sp_cmd->AddFlag( m_FastaFlag);
      // input for multiple files
      sp_cmd->AddFlag( m_FastaListFlag);
      // path for files
      sp_cmd->AddFlag( m_PathFlag);
      // flag for using pdb path hierarchy for all input files
      sp_cmd->AddFlag( m_PdbPathHierarchy);
      // path for output files
      sp_cmd->AddFlag( m_OutputPathFlag);
      // flag for calculating statistics
      sp_cmd->AddFlag( m_StatisticsFlag);
      // flag for reading predictions from different SSMethods
      sp_cmd->AddFlag( m_PdbListFlag);
      // flag for providing a oligomeric state dictionary
      sp_cmd->AddFlag( m_OligomericStateDictionaryFlag);
      // flag for adjusting the membrane
      sp_cmd->AddFlag( biol::Membrane::GetFlagMembrane());
      // flag indicating to define the path where membrane related data resides in not in the same path as m_PathFlag
      sp_cmd->AddFlag( m_MembraneDataPathFlag);

      // flag to convert to natural aa type
      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief default constructor
    SSPred::SSPred() :
      m_FastaFlag( new command::FlagStatic( "fasta", "flag for providing a single fasta sequence")),
      m_FastaParam( new command::Parameter( "fasta_tag", "id of sequence of interest, should not include any extensions", "")),
      m_FastaOligomericStateParam
      (
        new command::Parameter( "oligo", "oligomeric state should be 0 (monomer) or 1(multimer)", command::ParameterCheckRanged< size_t>( 0, 1), "0")
      ),
      m_FastaListFlag( new command::FlagStatic( "fasta_list", "flag for predictions for multiple sequence by providing a file that contains the list of fasta tags")),
      m_FastaListParam( new command::Parameter( "fasta_list_filename", "file that contains the list of fasta tags", ".")),
      m_PathFlag( new command::FlagStatic( "path", "flag for setting path where input files can be found")),
      m_PathParam( new command::Parameter( "path_param", "path where input files can be found", ".")),
      m_PdbPathHierarchy
      (
        new command::FlagStatic
        (
          "pdb_hierarchy",
          "boolean to indicate whether a pdb hiearchy is used in input paths so for pdbtag 1ABC.pdb it looks at {path}/AB/1ABC.pdb"
        )
      ),
      m_OutputPathFlag( new command::FlagStatic( "output_path", "flag for setting path where the output files should be written to")),
      m_OutputPathParam( new command::Parameter( "output_path_param", "path where the output files should be written to", ".")),
      m_StatisticsFlag
      (
        new command::FlagDynamic
        (
          "statistics", "flag for calculating statistics for the given SSMethods",
          command::Parameter
          (
            "ss_methods", "SSMethods of interest",
            command::ParameterCheckEnumerate< sspred::Methods>(),
            sspred::GetMethods().e_JUFO9D.GetName()
          ), 0, sspred::GetMethods().GetEnumCount()
        )
      ),
      m_PdbListFlag( new command::FlagStatic( "pdb_list", "flag to provide tags of pdb files to be read")),
      m_PdbListParam( new command::Parameter( "pdb_list_filename", "name of the file that has list of tags of pdbs to be read", "pdbs.ls")),
      m_OligomericStateDictionaryFlag( new command::FlagStatic( "oligo_dict", "flag for providing a dictionary for oligomeric states( pdbtag 0/1")),
      m_OligomericStateDictionaryParam( new command::Parameter( "oligo_dict_param", "file that contains the dictionary for oligomeric states", "")),
      m_MembraneDataPathFlag
      (
        new command::FlagStatic
        (
          "membrane_data_path",
          "flag indicating to define the path where membrane related data resides in not in the same path as m_PathFlag."
        )
      ),
      m_MembraneDataPathParam( new command::Parameter( "membrane_data_path_param", "path for where membrane related files should be searched for", ".")),
      m_PdbtmXmlFlag
      (
        new command::FlagStatic
        (
          "pdbtm_xml",
          "PDBTmXML files read from this path will be used to define a transformation that is applied to the pdb and contains the membrane core thickness,\
          that overwrites the core thickness given in the commandline "
        )
      )
    {
      // attach parameters to flags
      m_FastaFlag->PushBack( m_FastaParam);
      m_FastaFlag->PushBack( m_FastaOligomericStateParam);
      m_FastaListFlag->PushBack( m_FastaListParam);
      m_PathFlag->PushBack( m_PathParam);
      m_OutputPathFlag->PushBack( m_OutputPathParam);
      m_PdbListFlag->PushBack( m_PdbListParam);
      m_OligomericStateDictionaryFlag->PushBack( m_OligomericStateDictionaryParam);
      m_MembraneDataPathFlag->PushBack( m_MembraneDataPathParam);
    }

    const ApplicationType SSPred::SSPred_Instance
    (
      GetAppGroups().AddAppToGroup( new SSPred(), GetAppGroups().e_Sequence)
    );

  } // namespace app
} // namespace bcl
