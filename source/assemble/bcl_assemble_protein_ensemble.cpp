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
#include "assemble/bcl_assemble_protein_ensemble.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinEnsemble::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinEnsemble())
    );

    //! command line flag for reading list of protein models
    util::ShPtr< command::FlagInterface> &ProteinEnsemble::GetFlagPDBList()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "pdb_list", "\tlist of pdbs to do statistics over"
        )
      );

      // initialize filename parameter
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter
        (
          "filename", "\tfull path and name of pdb list file", "pdbs.ls"
        )
      );

      // initialize pdb column parameter
      static util::ShPtr< command::ParameterInterface> s_column
      (
        new command::Parameter
        (
          "column", "\tthe column in the file where the pdb name is. First column is number 0.", "0"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
        flag->PushBack( s_column);
      }

      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinEnsemble::ProteinEnsemble() :
      m_Ensemble()
    {
    }

    //! @brief constructor
    //! @param MODEL model the ensemble will be created from
    ProteinEnsemble::ProteinEnsemble( const util::ShPtr< ProteinModel> &MODEL) :
      ProteinModel( *util::ShPtr< ProteinModel>( MODEL->HardCopy())),
      m_Ensemble()
    {
    }

    //! @brief constructor
    //! @param MODEL model the ensemble will be created from
    ProteinEnsemble::ProteinEnsemble( const ProteinModel &MODEL) :
      ProteinModel( *util::ShPtr< ProteinModel>( MODEL.HardCopy())),
      m_Ensemble()
    {
    }

    //! @brief constructor taking filename with list of pdbs in it
    //! @param FILENAME file with the list of pdbs in it
    //! @param COLUMN the column in the file that has the pdbs in it
    //! @param AA_CLASS aa class used for protein models
    //! @param PREFIX prefix to add to all filenames in the file
    ProteinEnsemble::ProteinEnsemble
    (
      const std::string &FILENAME,
      const size_t COLUMN,
      const biol::AAClass &AA_CLASS,
      const std::string &PREFIX,
      const bool &STATUS_MESSAGES,
      const size_t &INPUT_START,
      const size_t &INPUT_MAX
    ) :
      m_Ensemble()
    {
      *this = GetEnsembleFromFile( FILENAME, COLUMN, AA_CLASS, PREFIX, STATUS_MESSAGES, INPUT_START, INPUT_MAX);
    }

    //! @brief Clone function
    //! @return pointer to new ProteinEnsemble
    ProteinEnsemble *ProteinEnsemble::Clone() const
    {
      return new ProteinEnsemble( *this);
    }

    //! @brief hard copy constructor
    //! @return a ProteinModel with chains hard copied from that model
    ProteinEnsemble *ProteinEnsemble::HardCopy() const
    {
      ProteinEnsemble *new_models( new ProteinEnsemble());
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        util::ShPtr< ProteinModel> new_model( ( *itr)->HardCopy());
        new_models->InsertElement( new_model);
        new_model->SetProteinModelData( ( *itr)->GetProteinModelData()->HardCopy());
      }

      // iterate over all chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( m_Chains.Begin()), chain_itr_end( m_Chains.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // push back a shptr to a hardcopied chain
        new_models->m_Chains.PushBack( util::ShPtr< Chain>( ( *chain_itr)->HardCopy()));
      }

      new_models->SetProteinModelData( m_ProteinModelData->HardCopy());

      return new_models;
    }

    //! @brief empty copy constructor
    //! @return a ProteinModel that is empty
    ProteinEnsemble *ProteinEnsemble::Empty() const
    {
      return new ProteinEnsemble();
    }

    //! @brief destructor
    ProteinEnsemble::~ProteinEnsemble()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns size of the container
    //! @return size, i.e. number of elements stored
    size_t ProteinEnsemble::GetSize() const
    {
      return m_Ensemble.GetSize();
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    ProteinEnsemble::iterator ProteinEnsemble::Begin()
    {
      return m_Ensemble.Begin();
    }

    //! @brief return const_iterator on begin
    //! @return const_iterator pointing to the beginning of the container, i.e. the first element
    ProteinEnsemble::const_iterator ProteinEnsemble::Begin() const
    {
      return m_Ensemble.Begin();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    ProteinEnsemble::iterator ProteinEnsemble::End()
    {
      return m_Ensemble.End();
    }

    //! @brief return const_iterator on end
    //! @return const_iterator pointing to the end of the container, i.e. behind the last element
    ProteinEnsemble::const_iterator ProteinEnsemble::End() const
    {
      return m_Ensemble.End();
    }

    //! @brief sets the sse pool within the protein model data
    //! @param SSE_POOL the sse pool that will be set in the protein model data
    void ProteinEnsemble::SetSSEPoolData( const util::ShPtr< SSEPool> &SSE_POOL)
    {
      // set the pool for the protein model
      ProteinModel::SetSSEPoolData( SSE_POOL);

      // iterate through the ensemble to set the pool
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        ( *itr)->SetSSEPoolData( SSE_POOL);
      }
    }

    //! @brief set identifiers of each model in the ensemble
    //! @param IDENTIFIERS the identifiers in the order of the models that they are assigned
    //! @return true if setting was successful
    bool ProteinEnsemble::SetIdentifiers( const storage::Vector< std::string> &IDENTIFIERS)
    {
      if( IDENTIFIERS.GetSize() != m_Ensemble.GetSize())
      {
        return false;
      }

      // iterate through the ensemble to set the pool
      storage::Vector< std::string>::const_iterator ident_itr( IDENTIFIERS.Begin());
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr, ++ident_itr)
      {
        util::ShPtr< ProteinModelData> sp_data( ( *itr)->GetProteinModelData());
        util::ShPtr< util::Wrapper< std::string> > sp_id( new util::Wrapper< std::string>( *ident_itr));
        if( !sp_data->Insert( ProteinModelData::e_Identification, sp_id))
        {
          return false;
        }
      }

      return true;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief checks whether container is empty
    //! @return if the container is empty or not
    bool ProteinEnsemble::IsEmpty() const
    {
      return m_Ensemble.IsEmpty();
    }

    //! @brief insert ELEMENT into the container
    //! @param ELEMENT an object of t_DataType that is inserted
    void ProteinEnsemble::InsertElement( const util::ShPtr< ProteinModel> &ELEMENT)
    {
      m_Ensemble.PushBack( ELEMENT);
    }

    //! @brief insert ELEMENT into the container
    //! @param POS the position where the element should be inserted
    //! @param ELEMENT an object of t_DataType that is inserted
    void ProteinEnsemble::InsertElement( const size_t POS, const util::ShPtr< ProteinModel> &ELEMENT)
    {
      m_Ensemble.InsertElement( POS, ELEMENT);
    }

    //! @brief delete single element at ITR
    //! @param ITR iterator pointing to the element that will be destroyed
    void ProteinEnsemble::RemoveElement( iterator ITR)
    {
      m_Ensemble.RemoveElement( ITR);
    }

    //! @brief gives the mean and std dev of distances for a data pair for this ensemble
    //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
    //! @return the mean and std dev of all pairwise distances
    math::RunningAverageSD< double> ProteinEnsemble::GetDistanceStatistics
    (
      const restraint::DataPairwise &DATA_PAIR
    ) const
    {
      // get the distances
      const storage::Vector< double> distances( GetDistances( DATA_PAIR));

      // will hold the statistics
      math::RunningAverageSD< double> mean_sd;

      // iterate through distances to calculate statistics
      for
      (
        storage::Vector< double>::const_iterator
          distances_itr( distances.Begin()), distances_itr_end( distances.End());
        distances_itr != distances_itr_end;
        ++distances_itr
      )
      {
        if( util::IsDefined( *distances_itr))
        {
          mean_sd += *distances_itr;
        }
      }

      // return statistics object
      return mean_sd;
    }

    //! @brief gives vector of coordinates located from each of the models in the ensemble
    //! @param COORD_LOCATOR method for locating coordinates in each model
    //! @return vector of coordinates located from each of the models in the ensemble
    storage::Vector< linal::Vector3D> ProteinEnsemble::GetCoordinates
    (
      const find::LocatorInterface< linal::Vector3D, ProteinModel> &COORD_LOCATOR
    ) const
    {
      // will hold coordinates
      storage::Vector< linal::Vector3D> coordinates;

      // iterate through the ensemble
      for
      (
        const_iterator model_itr( Begin()), model_itr_end( End()); model_itr != model_itr_end; ++model_itr
      )
      {
        // locate the coordinates
        const linal::Vector3D coord( COORD_LOCATOR.Locate( **model_itr));

        // true if the coordinates are not defined
        if( !coord.IsDefined())
        {
          // go to next iteration
          continue;
        }

        // add the coordinates to the vector of coordinates
        coordinates.PushBack( coord);
      }

      return coordinates;
    }

    //! @brief gives a vector of distances calculated within each of the models for a given data pair
    //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
    //! @return vector of doubles which are the distances calculated in each of the models of the ensemble
    storage::Vector< double> ProteinEnsemble::GetDistances( const restraint::DataPairwise &DATA_PAIR) const
    {
      // will hold the distances
      storage::Vector< double> distances;

      // iterate over the models of the ensemble
      for
      (
        const_iterator model_itr( Begin()), model_itr_end( End()); model_itr != model_itr_end; ++model_itr
      )
      {
        // get the coordinates for the first atom
        const linal::Vector3D coord_a( DATA_PAIR.First()->Locate( **model_itr));

        // get the coordinates for the second atom
        const linal::Vector3D coord_b( DATA_PAIR.Second()->Locate( **model_itr));

        // true if either of the coordinates are not defined
        if( !coord_a.IsDefined() || !coord_b.IsDefined())
        {
          // add undefined double and continue
          distances.PushBack( util::GetUndefinedDouble());

          continue;
        }

        // calculate the distance between the atoms
        const double distance( linal::Distance( coord_a, coord_b));

        // add the distance into distances
        distances.PushBack( distance);
      }

      return distances;
    }

    //! @brief gives a vector of distance changes for a data pair for this ensemble versus a given ensemble.
    //!        all pairwise distance changes are calculated and provided
    //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
    //! @param OTHER_ENSEMBLE the ensemble which provides the second state to get distance changes from
    //! @return vector of doubles which are all pairwise distance changes going from this ensemble to OTHER_ENSEMBLE
    //!         for the given data pair
    storage::Vector< double> ProteinEnsemble::GetDistanceChanges
    (
      const restraint::DataPairwise &DATA_PAIR,
      const ProteinEnsemble &OTHER_ENSEMBLE
    ) const
    {
      // get the distances from this ensemble
      const storage::Vector< double> this_distances( GetDistances( DATA_PAIR));

      // get the distances from OTHER_ENSEMBLE
      const storage::Vector< double> other_distances( OTHER_ENSEMBLE.GetDistances( DATA_PAIR));

      // vector will hold the pairwise distance changes between this and the other ensemble
      storage::Vector< double> distance_changes;

      // iterate through this distances to calculate the difference to the other distances
      for
      (
        storage::Vector< double>::const_iterator
          this_distances_itr( this_distances.Begin()), this_distances_itr_end( this_distances.End());
        this_distances_itr != this_distances_itr_end;
        ++this_distances_itr
      )
      {
        // iterate through the other distances
        for
        (
          storage::Vector< double>::const_iterator
            other_distances_itr( other_distances.Begin()), other_distances_itr_end( other_distances.End());
          other_distances_itr != other_distances_itr_end;
          ++other_distances_itr
        )
        {
          // calculate the distance difference
          const double distance_difference( ( *this_distances_itr) - ( *other_distances_itr));

          // add the distance difference to distance changes
          distance_changes.PushBack( distance_difference);
        }
      }

      return distance_changes;
    }

    //! @brief gives the mean and std dev of distance changes for a data pair for this ensemble versus a given
    //!        ensemble. all pairwise distance changes are calculated and provided
    //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
    //! @param OTHER_ENSEMBLE the ensemble which provides the second state to get distance changes from
    //! @return the mean and std dev of all pairwise distance changes going from this ensemble to OTHER_ENSEMBLE
    math::RunningAverageSD< double> ProteinEnsemble::GetDistanceChangesMeanSD
    (
      const restraint::DataPairwise &DATA_PAIR,
      const ProteinEnsemble &OTHER_ENSEMBLE
    ) const
    {
      // get the distance changes
      const storage::Vector< double> distance_changes( GetDistanceChanges( DATA_PAIR, OTHER_ENSEMBLE));

      // will hold the statistics
      math::RunningAverageSD< double> mean_sd;

      // iterate through distance changes to build up the mean and stdev statistics
      for
      (
        storage::Vector< double>::const_iterator
          distance_change_itr( distance_changes.Begin()), distance_change_itr_end( distance_changes.End());
        distance_change_itr != distance_change_itr_end;
        ++distance_change_itr
      )
      {
        // if the current distance change is not defined, continue
        if( !util::IsDefined( *distance_change_itr))
        {
          continue;
        }

        // add current value to statistics
        mean_sd += *distance_change_itr;
      }

      return mean_sd;
    }

    //! @brief gives all mean and std dev of distance changes for a data set for this ensemble versus a given
    //!        ensemble. all pairwise distance changes are calculated and provided
    //! @param DATA_SET the data set that indicates the distance that should be calculated
    //! @param OTHER_ENSEMBLE the ensemble which provides the second state to get distance changes from
    //! @return the mean and std dev of all pairwise distance changes going from this ensemble to OTHER_ENSEMBLE
    std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>
    ProteinEnsemble::GetDistanceChangesMeanSD
    (
      const restraint::DataSetPairwise &DATA_SET, const ProteinEnsemble &OTHER_ENSEMBLE
    ) const
    {
      // will hold the statistics for the whole data set
      std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean> stats;

      // iterate over the data set
      for
      (
        restraint::DataSetPairwise::const_iterator
          itr( DATA_SET.Begin()), itr_end( DATA_SET.End());
        itr != itr_end;
        ++itr
      )
      {
        // get and insert the current mean sd
        const math::RunningAverageSD< double> mean_sd( GetDistanceChangesMeanSD( *itr, OTHER_ENSEMBLE));
        stats.insert( std::make_pair( mean_sd, *itr));
      }

      return stats;
    }

    //! @brief removes all models from the ensemble
    void ProteinEnsemble::Reset()
    {
      m_Ensemble.Reset();
    }

    //! @brief read given ssmethods predictions
    //! @param SS_METHODS the methods to read for each protein in the ensemble
    //! @return true if reading was successful
    bool ProteinEnsemble::ReadSSPredictions( const storage::Set< sspred::Method> &SS_METHODS)
    {
      bool success( true);

      // iterate trough the ensemble
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        // get the filename
        util::ShPtr< util::Wrapper< std::string> > sp_filename( ( *itr)->GetProteinModelData()->GetData( ProteinModelData::e_PDBFile));

        if( !sp_filename.IsDefined())
        {
          success = false;
          continue;
        }

        // remove extension and split into path and prefix
        const std::string filename( io::File::RemoveFullExtension( *sp_filename));
        const storage::VectorND< 2, std::string> path_prefix( io::File::SplitToPathAndFileName( filename));

        // read predictions for that protein model
        success &=
          sspred::MethodHandler::ReadPredictionsForProteinModel
          (
            SS_METHODS,
            **itr,
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
          )
          || sspred::MethodHandler::ReadPredictionsForProteinModel( SS_METHODS, **itr, path_prefix.Second(), path_prefix.First());
      }

      // end
      return success;
    }

    //! @brief gives all of the pdb file names of the proteins in the ensemble
    //! @return vector of strings which are the pdb names of the proteins in the ensemble
    storage::Vector< std::string> ProteinEnsemble::GetPDBNames() const
    {
      storage::Vector< std::string> names;

      // iterate through the ensemble
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        // get the name of the first model in order to determine its length
        util::ShPtr< util::Wrapper< std::string> > pdb_name
        (
          ( *itr)->GetProteinModelData()->GetData( ProteinModelData::e_PDBFile)
        );
        if( pdb_name.IsDefined())
        {
          names.PushBack( std::string( *pdb_name));
        }
      }

      return names;
    }

    //! @brief gives a format object which can be used to format strings based on the longest pdb name
    //! @return format object which can be used to format strings to the same length as the longest pdb name
    util::Format ProteinEnsemble::GetNameFormatter() const
    {
      int longest_name( 0);
      // iterate through the ensemble
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        // get the name of the first model in order to determine its length
        util::ShPtr< util::Wrapper< std::string> > pdb_name
        (
          ( *itr)->GetProteinModelData()->GetData( ProteinModelData::e_PDBFile)
        );
        if( pdb_name.IsDefined())
        {
          int current_name_length( pdb_name->length());
          if( longest_name < current_name_length)
          {
            longest_name = current_name_length;
          }
        }
      }

      return util::Format().W( longest_name);
    }

    //! @brief gives the conformations that are in this protein
    //! @return protein ensemble which is the conformations making up this protein
    const ProteinEnsemble &ProteinEnsemble::GetConformationalEnsemble() const
    {
      return *this;
    }

    //! @brief sets the conformations that are in this protein
    //! @param ENSEMBLE protein ensemble which is the conformations making up this protein
    void ProteinEnsemble::SetConformationalEnsemble( const ProteinEnsemble &ENSEMBLE)
    {
      m_Ensemble = ENSEMBLE.m_Ensemble;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinEnsemble::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Ensemble, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Ensemble, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief opens a file and reads a list of pdbs and creates an ensemble out of them
    //! @param FILENAME file with the list of pdbs in it
    //! @param COLUMN the column in the file that has the pdbs in it - start at 0
    //! @param PREFIX prefix to add to all filenames in the file
    ProteinEnsemble ProteinEnsemble::GetEnsembleFromFile
    (
      const std::string &FILENAME,
      const size_t COLUMN,
      const biol::AAClass &AA_CLASS,
      const std::string &PREFIX,
      const bool &STATUS_MESSAGE,
      const size_t &INPUT_START,
      const size_t &INPUT_MAX
    )
    {
      // open file
      io::IFStream read;
      io::File::MustOpenIFStream( read, FILENAME);

      // read in and split the lines
      const storage::Vector< storage::Vector< std::string> > lines( util::SplittedStringLineListFromIStream( read));
      io::File::CloseClearFStream( read);
      ProteinEnsemble ensemble;

      const pdb::Factory factory( AA_CLASS);

      const size_t n_to_read( std::min( INPUT_MAX, lines.GetSize() - std::min( lines.GetSize(), INPUT_START)));

      const std::string n_pdbs( util::Format()( n_to_read));
      size_t pdb_num( 0);

      // iterate through the split lines to create protein models and add them to the ensemble
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          line_itr( lines.Begin() + INPUT_START), line_itr_end( lines.Begin() + n_to_read + INPUT_START);
        line_itr != line_itr_end;
        ++line_itr, ++pdb_num
      )
      {
        if( STATUS_MESSAGE)
        {
          util::GetLogger().LogStatus( "Reading Protein model : " + util::Format()( pdb_num) + " / " + n_pdbs);
        }
        // get pdb filename
        const std::string pdb_file( PREFIX + line_itr->operator()( COLUMN));

        // create protein model
        util::ShPtr< ProteinModel> model
        (
          factory.ProteinModelFromPDBFilename( pdb_file).Clone()
        );
        sspred::PDB::SetEnvironmentTypes( *model, true);

        // insert the pdb name as protein model data
        util::ShPtr< ProteinModelData> model_data( new ProteinModelData());
        model_data->Insert
        (
          ProteinModelData::e_PDBFile,
          util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( pdb_file))
        );
        model->SetProteinModelData( model_data);

        // add the model to the ensemble
        ensemble.InsertElement( model);
        BCL_MessageDbg( "inserted model " + pdb_file + " into ensemble");
      }

      return ensemble;
    }

  } // namespace assemble
} // namespace bcl
