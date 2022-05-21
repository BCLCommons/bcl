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
#include "assemble/bcl_assemble_sheet_template_handler.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_sheet.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_sse.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // static map to store the templates
    storage::Map< size_t, util::ShPtrVector< FoldTemplate> > SheetTemplateHandler::s_TemplateMap;

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SheetTemplateHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new SheetTemplateHandler())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new SheetTemplateHandler
    SheetTemplateHandler *SheetTemplateHandler::Clone() const
    {
      return new SheetTemplateHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SheetTemplateHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return command line flag for using sheet templates
    //! @return command line flag for using sheet templates
    util::ShPtr< command::FlagInterface> &SheetTemplateHandler::GetFlagSheetTemplates()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic( "sheet_templates", "\tuse sheet templates")
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_template_type
      (
        new command::Parameter
        (
          "sheet_template_type",
          "\ttype of sheet templates to use",
          command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "soluble", "membrane", "all")),
          "soluble"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_excluded_sheet_templates
      (
        new command::Parameter
        (
          "excluded_pdb_ids_file",
          "\tspecify a file that contains a list of PBD IDs to be excluded from the sheet template selection",
          ""
        )
      );

      // if this flag is initialized from the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_template_type);
        flag->PushBack( s_excluded_sheet_templates);
      }

      // end
      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief picks a random template of the appropriate size
    //! @param NR_STRANDS number of strands to be in the template
    //! @return a random template of the appropriate size
    const FoldTemplate &SheetTemplateHandler::GetRandomTemplate( const size_t NR_STRANDS)
    {
      // create static undefined template
      static const FoldTemplate s_undefined_template;

      // create an iterator to the vector of templates with the correct number of helices and strands
      storage::Map< size_t, util::ShPtrVector< FoldTemplate> >::iterator templates_itr
      (
        GetTemplates().Find( NR_STRANDS)
      );

      // if iterator is invalid, return empty fold template
      if( templates_itr == GetTemplates().End())
      {
        BCL_MessageStd
        (
          "Unable to find a sheet template with " + util::Format()( NR_STRANDS) + " strands"
        );
        return s_undefined_template;
      }

      // return a random template
      util::ShPtrVector< FoldTemplate>::iterator random_itr
      (
        random::GetGlobalRandom().Iterator
        (
          templates_itr->second.Begin(), templates_itr->second.End(), templates_itr->second.GetSize()
        )
      );

      BCL_MessageStd
      (
        "The sheet template is picked from " + ( *random_itr)->GetPDBID()
      )

      // if this is the first time this template is being picked
      // then the topology won't be initialized yet
      // therefore we need to collect the sheet topology and update the topology with this information
      if( !( *random_itr)->IsTopologyInitialized())
      {
        InitializeTopology( **random_itr);
      }

      // end
      return **random_itr;
    }

    //! @brief picks a random template of the appropriate size and geometry lengths
    //! @param SSES sses used to chose subtemplate based on length
    //! @param SSE_GEOMETRY_COMPARE comparison method
    //! @return a random template
    const FoldTemplate &SheetTemplateHandler::GetRandomTemplate
    (
      const util::SiPtrVector< const SSE> &SSES,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE
    )
    {
      // create static undefined template
      static const FoldTemplate s_undefined_template;

      // create an iterator to the vector of templates with the correct number of helices and strands
      storage::Map< size_t, util::ShPtrVector< FoldTemplate> >::iterator templates_itr( GetTemplates().Find( SSES.GetSize()));

      // if iterator is invalid, return empty fold template
      if( templates_itr == GetTemplates().End())
      {
        BCL_MessageStd
        (
          "Unable to find a sheet template with " + util::Format()( SSES.GetSize()) + " strands"
        );
        return s_undefined_template;
      }

      // initialize vector of matching templates
      util::ShPtrVector< FoldTemplate> matching_templates;

      // iterate through the templates
      for
      (
        util::ShPtrVector< FoldTemplate>::const_iterator template_itr( templates_itr->second.Begin()),
          template_itr_end( templates_itr->second.End());
        template_itr != template_itr_end; ++template_itr
      )
      {
        // if the geometries match the sses
        if( ( *template_itr)->HasSimilarSizeGeometries( SSES, SSE_GEOMETRY_COMPARE))
        {
          // add this template
          matching_templates.PushBack( *template_itr);
        }
      }

      // get a random template
      util::ShPtrVector< FoldTemplate>::iterator random_itr
      (
        random::GetGlobalRandom().Iterator
        (
          matching_templates.Begin(), matching_templates.End(), matching_templates.GetSize()
        )
      );

      // if no template was found
      if( random_itr == matching_templates.End())
      {
        BCL_MessageStd
        (
          "Unable to find a sheet template with " +
          util::Format()( SSES.GetSize()) + " strands that had the proper geometry sizes"
        );
        // return an empty template
        return s_undefined_template;
      }

      // if this is the first time this template is being picked
      // then the topology won't be initialized yet
      // therefore we need to collect the sheet topology and update the topology with this information
      if( !( *random_itr)->IsTopologyInitialized())
      {
        InitializeTopology( **random_itr);
      }

      // return a random template
      return **random_itr;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SheetTemplateHandler::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SheetTemplateHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief singleton function to return sheet templates previously generated and stored in a file
    //! @return sheet templates previously generated and stored in a file
    storage::Map< size_t, util::ShPtrVector< FoldTemplate> > &SheetTemplateHandler::GetTemplates()
    {
      // if the map is not initialized
      if( s_TemplateMap.IsEmpty())
      {
        // get the template type
        const std::string &template_type( GetFlagSheetTemplates()->GetFirstParameter()->GetValue());

        // if both soluble and membrane are to be used
        if( template_type == "all")
        {
          // read in both
          s_TemplateMap = ReadTemplates
          (
            storage::Set< std::string>( GetInputFilename( "soluble"), GetInputFilename( "membrane"))
          );
        }
        // only one type is to be used
        else
        {
          // read in that type
          s_TemplateMap = ReadTemplates( storage::Set< std::string>( GetInputFilename( template_type)));
        }
      }

      //end
      return s_TemplateMap;
    }

    //! @brief singleton function to return sheet templates previously generated and stored in a file
    //! @param FILENAMES set of filenames to be read in
    //! @return sheet templates previously generated and stored in a file
    storage::Map< size_t, util::ShPtrVector< FoldTemplate> > SheetTemplateHandler::ReadTemplates
    (
      const storage::Set< std::string> &FILENAMES
    )
    {
      // initialize storage map
      storage::Map< size_t, util::ShPtrVector< FoldTemplate> > sheet_template_map;

      // initialize collector
      const util::ShPtr< CollectorTopologyInterface> sp_collector
      (
        new CollectorTopologySheet()
      );

      // iterate through the filenames
      for
      (
        storage::Set< std::string>::const_iterator file_itr( FILENAMES.Begin()), file_itr_end( FILENAMES.End());
        file_itr != file_itr_end; ++file_itr
      )
      {
        // open input file
        io::IFStream read;
        io::File::MustOpenIFStream( read, *file_itr);

        // read in the number of templates
        size_t number_of_templates;
        read >> number_of_templates;

        // iterate through the number of templates
        for( size_t i( 0); i != number_of_templates; ++i)
        {
          // initialize vectors of geometries
          util::ShPtrVector< SSEGeometryPhiPsi> geometry_vector;

          // read the number of geometries and the pdb id
          size_t nr_strands;
          std::string pdb_id;
          read >> pdb_id >> nr_strands;

          // iterate over the number of geometries
          for( size_t j( 0); j != nr_strands; ++j)
          {
            // read in the type
            std::string ss_type_string;
            read >> ss_type_string;
            biol::SSType ss_type( biol::GetSSTypes().GetEnumFromName( ss_type_string));

            // read in N coords
            linal::Vector3D n_coords;
            for( size_t k( 0); k != 3; ++k)
            {
              read >> n_coords( k);
            }

            // read in CA coords
            linal::Vector3D ca_coords;
            for( size_t k( 0); k != 3; ++k)
            {
              read >> ca_coords( k);
            }

            // read in C coords
            linal::Vector3D c_coords;
            for( size_t k( 0); k != 3; ++k)
            {
              read >> c_coords( k);
            }

            // read in phi/psi size
            size_t nr_residues;
            read >> nr_residues;

            // read in the first psi
            storage::Vector< storage::VectorND< 2, double> > angles;
            double first_psi;
            read >> first_psi;
            angles.PushBack( storage::VectorND< 2, double>( util::GetUndefined< double>(), first_psi));

            // iterate over the number of residues ( minus 2 since the first and last angles are nan) to read in phi/psis
            for( size_t l( 0); l != nr_residues - 2; ++l)
            {
              double phi;
              double psi;
              read >> phi >> psi;

              angles.PushBack( storage::VectorND< 2, double>( phi, psi));
            }

            // read in the last phi
            double last_phi;
            read >> last_phi;
            angles.PushBack( storage::VectorND< 2, double>( last_phi, util::GetUndefined< double>()));

            // add the geometry to the vector
            geometry_vector.PushBack
            (
              util::ShPtr< SSEGeometryPhiPsi>
              (
                new SSEGeometryPhiPsi
                (
                  biol::AASequencePhiPsi( n_coords, ca_coords, c_coords, angles),
                  ss_type,
                  "S" + util::Format()( nr_strands)
                )
              )
            );
          }

          // if the template is not to be excluded
          if( !GetExcludedPdbs().Contains( pdb_id))
          {
            // create the fold template and add it to the map
            sheet_template_map[ nr_strands].PushBack
            (
              util::ShPtr< FoldTemplate>( new FoldTemplate( geometry_vector, sp_collector, pdb_id, false))
            );
          }
        }

        // close the read stream
        io::File::CloseClearFStream( read);
      }

      // end
      return sheet_template_map;
    }

    //! @brief default input filename
    //! @param TEMPLATE_TYPE type of templates to use
    //! @return default input filename
    const std::string &SheetTemplateHandler::GetInputFilename( const std::string &TEMPLATE_TYPE)
    {
      // initialize with default file names
      static const std::string s_soluble_filename( score::Score::AddHistogramPath( "soluble_sheet_templates.input"));
      static const std::string s_membrane_filename( score::Score::AddHistogramPath( "membrane_sheet_templates.input"));
      static const std::string s_undefined;

      // return filename determined by fold template type
      if( TEMPLATE_TYPE == "soluble")
      {
        return s_soluble_filename;
      }
      else if( TEMPLATE_TYPE == "membrane")
      {
        return s_membrane_filename;
      }

      // type should be soluble or membrane so if this point is reached return default string
      return s_undefined;
    }

    //! @brief return set of excluded pdb ids
    //! @return set of excluded pdb ids
    storage::Set< std::string> &SheetTemplateHandler::GetExcludedPdbs()
    {
      static storage::Set< std::string> s_excluded_pdbs;

      // if a file that contains PDBs to exclude as templates was passed
      if( s_excluded_pdbs.IsEmpty() && GetFlagSheetTemplates()->GetParameterList()( 1)->GetWasSetInCommandLine())
      {
        // open the file
        io::IFStream read;
        io::File::MustOpenIFStream( read, GetFlagSheetTemplates()->GetParameterList()( 1)->GetValue());

        // read in the pdb ids
        while( !read.eof())
        {
          std::string temp_string;
          read >> temp_string;

          // add the string to the set
          if( temp_string != "")
          {
            std::transform( temp_string.begin(), temp_string.end(), temp_string.begin(), toupper);
            s_excluded_pdbs.Insert( temp_string);
          }
        }
        io::File::CloseClearFStream( read);
      }

      // end
      return s_excluded_pdbs;
    }

    //! @brief initialize the sheet template topology
    //! @param SHEET_TEMPLATE sheet template for which topology will be initialized
    void SheetTemplateHandler::InitializeTopology( FoldTemplate &SHEET_TEMPLATE)
    {
      // calculate the geometries
      SHEET_TEMPLATE.CalculateSSEGeometries();

      // collect the sheet topology from the given vector of geometries
      util::ShPtrVector< Topology> sheet_vector( CollectorTopologySheet().Collect( SHEET_TEMPLATE.GetGeometries()));

      BCL_MessageStd
      (
        "The picked template, has " + util::Format()( sheet_vector.GetSize()) +
        " sheet/s and #strands for first one: " + util::Format()( sheet_vector.FirstElement()->GetElements().GetSize())
      );

      // make sure there is only one sheet and is of correct size
      if
      (
        sheet_vector.GetSize() != 1 ||
        sheet_vector.FirstElement()->GetElements().GetSize() != SHEET_TEMPLATE.GetGeometries().GetSize()
      )
      {
        BCL_MessageCrt( "There was a problem with sheet topology building from sheet template selected");
      }
      else
      {
        // update the sheet topology
        SHEET_TEMPLATE.SetTopology( *sheet_vector.FirstElement());
      }
    }

  } // namespace assemble
} // namespace bcl
