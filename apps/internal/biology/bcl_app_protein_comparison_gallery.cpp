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
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_wrapper.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinComparisonGallery
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_app_ProteinComparisonGallery.cpp @endlink
    //! @author alexanns, weinerbe
    //! @date Mar 8, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinComparisonGallery :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

        util::ShPtr< command::FlagInterface> m_SortFlag;          //!< sorting flag
        util::ShPtr< command::FlagInterface> m_SortCriteriaFlag;  //!< sorting method param
        util::ShPtr< command::FlagInterface> m_ScoreTableFlag;    //!< score table flag
        util::ShPtr< command::FlagInterface> m_FlagNumPicsToMake; //!< how many pictures of proteins to make

        util::ShPtr< command::FlagInterface> m_PixelsX; //!< size of the generated png
        util::ShPtr< command::FlagInterface> m_PixelsY; //!< size of the generated png
        util::ShPtr< command::FlagInterface> m_OutputFilename; //!< filename of the script to be generated

        //! flag indicating it is a membrane protein and should be rotated for better viewing
        util::ShPtr< command::FlagInterface> m_MembraneOrient;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinComparisonGallery();

      //! @brief Clone function
      //! @return pointer to new ProteinComparisonGallery
      ProteinComparisonGallery *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the Command object
      util::ShPtr< command::Command> InitializeCommand() const;

    ////////////////
    // operations //
    ////////////////

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType ProteinComparisonGallery_Instance;

    }; // class ProteinComparisonGallery

    //! @brief default constructor
    ProteinComparisonGallery::ProteinComparisonGallery() :
      m_SortFlag
      (
        new command::FlagDynamic
        (
          "sort",
          "sort score table with by column name(s), if multiple given, the sum will be considered",
          command::Parameter( "sort_column", "name of column to be sorted")
        )
      ),
      m_SortCriteriaFlag
      (
        new command::FlagStatic
        (
          "sort_criteria",
          "sorting criteria",
          command::Parameter
          (
            "sort_criteria",
            "sorting criteria",
            command::ParameterCheckEnumerate< math::Comparisons< double> >(),
            math::Comparisons< double>::GetEnums().e_Less.GetName()
          )
        )
      ),
      m_ScoreTableFlag
      (
        new command::FlagStatic
        (
          "table",
          "score table to be read in",
          command::Parameter( "table", "score table to be read in", "")
        )
      ),
      m_FlagNumPicsToMake
      (
        new command::FlagStatic
        (
          "num_pics",
          "the number of pictures that should be made. one per model of the best models accoring to sort criteria",
          command::Parameter( "number of pics", "how many models to see", "1")
        )
      ),
      m_PixelsX
      (
        new command::FlagStatic
        (
          "pixels_x",
          "resolution in the x-direction",
          command::Parameter( "horizontal resolution", "size_t which is horizontal resolution", "900")
        )
      ),
      m_PixelsY
      (
        new command::FlagStatic
        (
          "pixels_y",
          "resolution in the y-direction",
          command::Parameter( "vertical resolution", "size_t which is vertical resolution", "600")
        )
      ),
      m_OutputFilename
      (
        new command::FlagStatic
        (
          "output_filename",
          "the name of the outputted pymol script for generating the pngs",
          command::Parameter( "output_filename", "string which is the file the script will be written to", "make_pngs.pml")
        )
      ),
      m_MembraneOrient
      (
        new command::FlagStatic
        (
          "membrane_orient",
          "Indicates this is a membrane protein and it should be rotated to improve the view"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldAnalysis
      ProteinComparisonGallery *ProteinComparisonGallery::Clone() const
    {
      return new ProteinComparisonGallery( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinComparisonGallery::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the Command object
    util::ShPtr< command::Command> ProteinComparisonGallery::InitializeCommand() const
    {
      // create command
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member flags
      sp_cmd->AddFlag( m_SortFlag);
      sp_cmd->AddFlag( m_SortCriteriaFlag);
      sp_cmd->AddFlag( m_ScoreTableFlag);
      sp_cmd->AddFlag( m_FlagNumPicsToMake);
      sp_cmd->AddFlag( m_PixelsX);
      sp_cmd->AddFlag( m_PixelsY);
      sp_cmd->AddFlag( m_OutputFilename);
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagPrefix());
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagNativeModel());
      sp_cmd->AddFlag( m_MembraneOrient);

      // add storage flag
      sp_cmd->AddFlag( assemble::ProteinStorageFile::GetDefaultStorageFlag());

      // add default bcl flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled command
      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int ProteinComparisonGallery::Main() const
    {
      // get all the table filenames
      std::string table_filename( m_ScoreTableFlag->GetFirstParameter()->GetValue());

      // table
      storage::Table< double> table;

      // read it in
      io::IFStream read;
      io::File::MustOpenIFStream( read, table_filename);
      table.ReadFormatted( read);
      io::File::CloseClearFStream( read);

      // table with desired columns
      storage::Table< double> tables_columns;

      // copy the header
      storage::TableHeader new_header( table.GetHeader());

      // add a new column
      const std::string sum_column_name( "sum_for_sort");
      new_header.PushBack( sum_column_name);

      // table for sorted scores
      storage::Table< double> table_sorted_scores( new_header);

      // initialize vector of indices
      storage::Vector< size_t> indices;

      // iterate over the passed column names
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          param_itr( m_SortFlag->GetParameterList().Begin()),
          param_itr_end( m_SortFlag->GetParameterList().End());
        param_itr != param_itr_end; ++param_itr
      )
      {
        // check if the table has this column
        if( new_header.HasColumn( ( *param_itr)->GetValue()))
        {
          indices.PushBack( new_header[ ( *param_itr)->GetValue()]);
        }
        else
        {
          BCL_MessageStd
          (
            "No column found for " + ( *param_itr)->GetValue() + ", so skipping it for sorting"
          );
        }
      }

      // iterate over the original table
      for
      (
        storage::Table< double>::const_iterator table_itr( table.Begin()), table_itr_end( table.End());
        table_itr != table_itr_end; ++table_itr
      )
      {
        // get the data from the row
        storage::Vector< double> row_data( table_itr->Second().GetData());

        // initialize sum
        double sum( 0);

        // iterate over the indices
        for
        (
          storage::Vector< size_t>::const_iterator index_itr( indices.Begin()), index_itr_end( indices.End());
          index_itr != index_itr_end; ++index_itr
        )
        {
          sum += row_data( *index_itr);
        }

        // add the sum to the vector
        row_data.PushBack( sum);

        // insert a new row
        table_sorted_scores.InsertRow( table_itr->First(), row_data);
      }

      // now sort the table by the new sum column
      table_sorted_scores.SortByColumn
      (
        sum_column_name,
        **math::Comparisons< double>::Comparison( m_SortCriteriaFlag->GetFirstParameter()->GetValue())
      );

      assemble::ProteinEnsemble proteins;

      // get the storage object
      util::ShPtr< assemble::ProteinStorageFile> storage( assemble::ProteinStorageFile::GetDefaultStorage());
      BCL_Assert( storage.IsDefined(), "storage is not defined properly");

      const std::string output_prefix( fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue());

      // to hold values for the rows
      storage::Vector< double> sum_column_values;

      // iterate through the table to get the desired number of pdb names
      size_t counter( 0);
      for
      (
        storage::Table< double>::const_iterator
          table_itr( table_sorted_scores.Begin()), table_itr_end( table_sorted_scores.End());
        table_itr != table_itr_end && counter < m_FlagNumPicsToMake->GetFirstParameter()->GetNumericalValue< size_t>();
        ++table_itr, ++counter
      )
      {
        const double sum_column_value( table_itr->Second()[ sum_column_name]);
        sum_column_values.PushBack( sum_column_value);

        // protein storage specifies path
        if( assemble::ProteinStorageFile::GetDefaultStorageFlag()->GetFlag())
        {
          io::DirectoryEntry file( table_itr->First());
          BCL_MessageDbg( "full name " + file.GetFullName());
          std::string name( io::File::RemoveFullExtension( file.GetName()));
          BCL_MessageDbg( "protein storage will use name " + name);

          // create the name of the pdb that will be outputted
          std::string output_pdb_filename
          (
            output_prefix + name + "_" + util::Format()( sum_column_value) + ".pdb"
          );

          // check for key sizes = s_SmallKeySize (final models)
          if( name.size() > assemble::ProteinStorageFile::s_SmallKeySize)
          {
            // get the substring
            const std::string key_short
            (
              name.substr
              (
                name.length() - assemble::ProteinStorageFile::s_SmallKeySize,
                assemble::ProteinStorageFile::s_SmallKeySize
              )
            );

            // if this key is valid
            if( assemble::ProteinStorageFile::IsValidKey( key_short))
            {
              util::ShPtr< assemble::ProteinModel> model
              (
                storage->Retrieve
                (
                  name.substr( 0, name.length() - assemble::ProteinStorageFile::s_SmallKeySize), key_short
                )
              );
              BCL_Assert( model.IsDefined(), "could not retrieve model " + name);
              util::ShPtr< assemble::ProteinModelData> protein_model_data( model->GetProteinModelData());
              protein_model_data->Insert
              (
                assemble::ProteinModelData::e_PDBFile,
                util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( output_pdb_filename))
              );
              proteins.InsertElement( model);
            }

            // check for key sizes = s_LargeKeySize
            else if( name.size() > assemble::ProteinStorageFile::s_LargeKeySize)
            {
              // get the substring
              const std::string key_long
              (
                name.substr
                (
                  name.length() - assemble::ProteinStorageFile::s_LargeKeySize,
                  assemble::ProteinStorageFile::s_LargeKeySize
                )
              );

              // if this key is valid
              if( assemble::ProteinStorageFile::IsValidKey( key_long))
              {
                util::ShPtr< assemble::ProteinModel> model
                (
                  storage->Retrieve
                  (
                    name.substr( 0, name.length() - assemble::ProteinStorageFile::s_LargeKeySize), key_long
                  )
                );
                BCL_Assert( model.IsDefined(), "could not retrieve model " + name);
                util::ShPtr< assemble::ProteinModelData> protein_model_data( model->GetProteinModelData());
                protein_model_data->Insert
                (
                  assemble::ProteinModelData::e_PDBFile,
                  util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( output_pdb_filename))
                );
                proteins.InsertElement( model);
              }
            }
          }
        }
        else //< protein storage not specified, assume full path and filename are the row name
        {
          std::string name( table_itr->First());
          BCL_MessageDbg( "full name " + name);
          pdb::Factory factory;
          util::ShPtr< assemble::ProteinModel> model
          (
            factory.ProteinModelFromPDBFilename( name).Clone()
          );
          BCL_Assert( model.IsDefined(), "could not retrieve model " + name);
          util::ShPtr< assemble::ProteinModelData> protein_model_data( model->GetProteinModelData());

          // create the name of the pdb that will be outputted
          io::DirectoryEntry file( name);
          std::string output_pdb_filename
          (
            output_prefix + io::File::RemoveFullExtension( io::File::RemoveFullExtension( file.GetName())) +
            "_" + util::Format()( sum_column_value) + ".pdb"
          );

          BCL_MessageDbg( "writing to file " + output_pdb_filename);

          protein_model_data->Insert
          (
            assemble::ProteinModelData::e_PDBFile,
            util::ShPtr< util::Wrapper< std::string> >( new util::Wrapper< std::string>( output_pdb_filename))
          );
          proteins.InsertElement( model);
        }
      }

      BCL_MessageDbg( "protein size is " + util::Format()( proteins.GetSize()));
      // get the native model
      // get the pdb tag from the row name
      const std::string native_pdb_filename( fold::DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue());

      // create the protein model
      pdb::Factory factory;
      const assemble::ProteinModel native_model( factory.ProteinModelFromPDBFilename( native_pdb_filename));

      // get the pymol output script filename
      std::string output_filename( m_OutputFilename->GetFirstParameter()->GetValue());

      // open the output pymol script
      io::OFStream script;
      io::File::MustOpenOFStream( script, output_filename);

      // load the native models
      script << "load " << native_pdb_filename << ", native\n";
      script << "set_color greytrans, [0.7, 0.7, 0.7]\n";
      script << "set depth_cue = 0\n";
      script << "hide everything, native\n";
      script << "show cartoon, native\n";
      script << "set cartoon_transparency, 0.3, native\n";
      script << "color greytrans, native\n";

      if( m_MembraneOrient->GetFlag())
      {
        script << "rotate x, 90\n";
      }

      script << "bg white\n";

      size_t model_counter( 0);
      pdb::Factory writing_factory;
      // iterate over the models and load them
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( proteins.Begin()), ensemble_itr_end( proteins.End());
          ensemble_itr != ensemble_itr_end;
        ++ensemble_itr, ++model_counter
      )
      {
        // get copy of current model
        const assemble::ProteinModel &current_model( **ensemble_itr);
        util::ShPtr< util::Wrapper< std::string> > name
        (
          current_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );
        BCL_Assert( name.IsDefined(), "could not get protein model data of pdb name from model");
        std::string pymol_name( "model_" + util::Format()( model_counter));
        script << "load " << name->GetData() << ", " << pymol_name << "\n";
        script << "hide everything," << pymol_name << "\n";
        script << "show cartoon," << pymol_name << "\n";
        script << "spectrum count, selection =" << pymol_name << ", byres=1\n";
        if( m_MembraneOrient->GetFlag())
        {
          script << "rotate x, 90, " << pymol_name << "\n";
        }
        script << "zoom native\n";
        script << "ray " << m_PixelsX->GetFirstParameter()->GetNumericalValue< size_t>()
               << ","    << m_PixelsY->GetFirstParameter()->GetNumericalValue< size_t>()
               << "\n";
        const std::string png_name( pymol_name + "_" + util::Format()( sum_column_values( model_counter)) + ".png");
        script << "png " << output_prefix + png_name << "\n";
        script << "hide everything, " << pymol_name << "\n";

        io::OFStream current_write;
        io::File::MustOpenOFStream( current_write, name->GetData());
        writing_factory.WriteModelToPDB( current_model, current_write);
        io::File::CloseClearFStream( current_write);
      }
      script << "quit\n";
      io::File::CloseClearFStream( script);

      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinComparisonGallery::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinComparisonGallery::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    const ApplicationType ProteinComparisonGallery::ProteinComparisonGallery_Instance( GetAppGroups().AddAppToGroup( new ProteinComparisonGallery(), GetAppGroups().e_Protein));

  } // namespace app
} // namespace bcl
