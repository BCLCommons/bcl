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
#include "assemble/bcl_assemble_fold_template_handler.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_sse.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // static map to store the templates
    storage::Map< storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate> > FoldTemplateHandler::s_TemplateMap;

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FoldTemplateHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new FoldTemplateHandler())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FoldTemplateHandler::FoldTemplateHandler()
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldTemplateHandler
    FoldTemplateHandler *FoldTemplateHandler::Clone() const
    {
      return new FoldTemplateHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FoldTemplateHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return command line flag for using fold templates
    //! @return command line flag for using fold templates
    util::ShPtr< command::FlagInterface> &FoldTemplateHandler::GetFlagFoldTemplates()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic( "fold_templates", "\tfit protein models to fold templates")
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_template_type
      (
        new command::Parameter
        (
          "fold_template_type",
          "\ttype of fold templates to use",
          command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "soluble", "membrane", "all")),
          "soluble"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_small_template_probability
      (
        new command::Parameter
        (
          "small_template_probability",
          "\tprobability of using a small template (less bodies than SSEs in the pool) for the initial placement",
          "1"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_standard_template_probability
      (
        new command::Parameter
        (
          "standard_template_probability",
          "\tprobability of using a standard template (same bodies as SSEs in the pool) for the initial placement",
          "1"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_large_template_probability
      (
        new command::Parameter
        (
          "large_template_probability",
          "\tprobability of using a large template (more bodies than SSEs in the pool) for the initial placement",
          "1"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_excluded_fold_templates
      (
        new command::Parameter
        (
          "excluded_pdb_ids_file",
          "\tspecify a file that contains a list of PBD IDs to be excluded from the fold template selection",
          ""
        )
      );

      // if this flag is initialized from the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_template_type);
        flag->PushBack( s_small_template_probability);
        flag->PushBack( s_standard_template_probability);
        flag->PushBack( s_large_template_probability);
        flag->PushBack( s_excluded_fold_templates);
      }

      // end
      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief picks a random template of the appropriate size
    //! @param HELICES number of helices to be in the template
    //! @param STRANDS number of strands to be in the template
    //! @return a random template of the appropriate size
    const FoldTemplate &FoldTemplateHandler::GetRandomTemplate( const size_t HELICES, const size_t STRANDS)
    {
      // create static undefined template
      static const FoldTemplate s_undefined_template;

      // create an iterator to the vector of templates with the correct number of helices and strands
      storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      >::const_iterator templates_itr
      (
        GetTemplates().Find( storage::Pair< size_t, size_t>( HELICES, STRANDS))
      );

      // if iterator is invalid, return empty fold template
      if( templates_itr == GetTemplates().End())
      {
        BCL_MessageStd
        (
          "Unable to find a fold template with " + util::Format()( HELICES) +
            " helices and " + util::Format()( STRANDS) + " strands"
        );
        return s_undefined_template;
      }

      // return a random template
      util::ShPtrVector< FoldTemplate>::const_iterator random_itr
      (
        random::GetGlobalRandom().Iterator
        (
          templates_itr->second.Begin(), templates_itr->second.End(), templates_itr->second.GetSize()
        )
      );

      // end
      return **random_itr;
    }

    //! @brief picks a random template of the appropriate size and geometry lengths
    //! @param SSES sses used to chose subtemplate based on length
    //! @param SSE_GEOMETRY_COMPARE comparison method
    //! @return a random template
    const FoldTemplate &FoldTemplateHandler::GetRandomTemplate
    (
      const util::SiPtrVector< const SSE> &SSES,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE
    )
    {
      // get the helical sses
      const CollectorSSE helix_collector( biol::GetSSTypes().HELIX);
      util::SiPtrList< const SSE> helices( helix_collector.Collect( SSES));

      // get the strand sses
      const CollectorSSE strand_collector( biol::GetSSTypes().STRAND);
      util::SiPtrList< const SSE> strands( strand_collector.Collect( SSES));

      // create static undefined template
      static const FoldTemplate s_undefined_template;

      // create an iterator to the vector of templates with the correct number of helices and strands
      storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      >::const_iterator templates_itr
      (
        GetTemplates().Find( storage::Pair< size_t, size_t>( helices.GetSize(), strands.GetSize()))
      );

      // if iterator is invalid, return empty fold template
      if( templates_itr == GetTemplates().End())
      {
        BCL_MessageStd
        (
          "Unable to find a fold template with " + util::Format()( helices.GetSize()) +
            " helices and " + util::Format()( strands.GetSize()) + " strands"
        );
        return s_undefined_template;
      }

      BCL_MessageDbg( "# potential templates: " + util::Format()( templates_itr->second.GetSize()));

      // initialize vector of matching templates
      util::SiPtrVector< const FoldTemplate> matching_templates;

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
      util::SiPtrVector< const FoldTemplate>::const_iterator random_itr
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
          "Unable to find a fold template with " + util::Format()( helices.GetSize()) +
            " helices and " + util::Format()( strands.GetSize()) + " strands that had the proper geometry sizes"
        );
        // return an empty template
        return s_undefined_template;
      }

      // return a random template
      return **random_itr;
    }

    //! @brief returns a random fold template (that has geometries close in size to the passed sses)
    //!        generated from a larger template
    //! @param SSES sses used to chose subtemplate based on length
    //! @param SSE_GEOMETRY_COMPARE comparison method
    //! @return a random fold template generated from a larger template
    FoldTemplate FoldTemplateHandler::GetRandomSubTemplate
    (
      const util::SiPtrVector< const SSE> &SSES,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE
    )
    {
      // get the helical sses
      const CollectorSSE helix_collector( biol::GetSSTypes().HELIX);
      util::SiPtrList< const SSE> helices( helix_collector.Collect( SSES));

      // get the strand sses
      const CollectorSSE strand_collector( biol::GetSSTypes().STRAND);
      util::SiPtrList< const SSE> strands( strand_collector.Collect( SSES));

      // get the large templates
      util::ShPtrVector< FoldTemplate> large_templates( GetLargeTemplates( helices, strands, SSE_GEOMETRY_COMPARE));

      // initialize sub-template
      FoldTemplate sub_template;

      // loop until a sub template is found
      while( sub_template.GetGeometries().IsEmpty() && !large_templates.IsEmpty())
      {
        // get a random template iterator
        util::ShPtrVector< FoldTemplate>::iterator large_template_itr
        (
          random::GetGlobalRandom().Iterator( large_templates.Begin(), large_templates.End(), large_templates.GetSize())
        );

        // if the topology has not been initialized
        if( !( *large_template_itr)->IsTopologyInitialized())
        {
          // calculate the topology
          ( *large_template_itr)->CalculateTopology();
        }

        // get a subdomain and set it to the sub_template
        sub_template = ( *large_template_itr)->GetSubDomain( SSES, SSE_GEOMETRY_COMPARE);

        // remove the template from the vector so it is not used again if no subdomains were found
        large_templates.RemoveElement( large_template_itr);
      }

      // the sub-template is still empty
      if( sub_template.GetGeometries().IsEmpty())
      {
        // warn the user
        BCL_MessageStd
        (
          "Unable to find a subtemplate, trying a full template"
        );

        // get a random full template
        GetRandomTemplate( helices.GetSize(), strands.GetSize());
      }

      // end
      return sub_template;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FoldTemplateHandler::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FoldTemplateHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief singleton function to return fold templates previously generated and stored in a file
    //! @return fold templates previously generated and stored in a file
    const storage::Map
    <
      storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
    > &FoldTemplateHandler::GetTemplates()
    {
      // if the map is not initialized
      if( s_TemplateMap.IsEmpty())
      {
        // get the template type
        const std::string &template_type( GetFlagFoldTemplates()->GetFirstParameter()->GetValue());

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

    //! @brief singleton function to return fold templates previously generated and stored in a file
    //! @param FILENAMES set of filenames to be read in
    //! @return fold templates previously generated and stored in a file
    storage::Map
    <
      storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
    > FoldTemplateHandler::ReadTemplates( const storage::Set< std::string> &FILENAMES)
    {
      // initialize storage map
      storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      > fold_template_map;

      // initialize collector
      const util::ShPtr< CollectorTopologyInterface> sp_collector
      (
        new CollectorTopologyCombined()
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
          size_t nr_geometries;
          std::string pdb_id;
          read >> pdb_id >> nr_geometries;

          geometry_vector.AllocateMemory( nr_geometries);

          // initialize helix count
          size_t helix_ctr( 0);
          size_t strand_ctr( 0);

          // iterate over the number of geometries
          for( size_t j( 0); j != nr_geometries; ++j)
          {
            // read in the type
            std::string ss_type_string;
            read >> ss_type_string;
            biol::SSType ss_type( biol::GetSSTypes().GetEnumFromName( ss_type_string));

            // initialize identification
            std::string identification;

            // if this is a helix
            if( ss_type == biol::GetSSTypes().HELIX)
            {
              identification = "H" + util::Format()( helix_ctr);
              ++helix_ctr;
            }
            else
            {
              identification = "S" + util::Format()( strand_ctr);
              ++strand_ctr;
            }

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
            angles.AllocateMemory( nr_residues);
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
                  identification
                )
              )
            );
          }

          // if the template is not to be excluded
          if( !GetExcludedPdbs().Contains( pdb_id))
          {
            // create the fold template and add it to the map
            fold_template_map[ storage::Pair< size_t, size_t>( helix_ctr, strand_ctr)].PushBack
            (
              util::ShPtr< FoldTemplate>( new FoldTemplate( geometry_vector, sp_collector, pdb_id))
            );
          }
          else
          {
            BCL_MessageStd( "Ignoring fold template " + pdb_id);
          }
        }

        // close the read stream
        io::File::CloseClearFStream( read);
      }

      // end
      return fold_template_map;
    }

    //! @brief gets all fold templates from the requested composition up to a constant size
    //! @param HELICES helical SSEs to be used to find a suitable template
    //! @param STRANDS strand SSEs to be used to find a suitable template
    //! @param SSE_GEOMETRY_COMPARE comparison method
    //! @return all fold templates from the requested composition up to a constant size
    util::ShPtrVector< FoldTemplate> FoldTemplateHandler::GetLargeTemplates
    (
      const util::SiPtrList< const SSE> &HELICES,
      const util::SiPtrList< const SSE> &STRANDS,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE
    )
    {
      // max number of helices or strands possible
      static const size_t s_max_sses( 10);

      // initialize vector of potential templates
      util::ShPtrVector< FoldTemplate> large_templates;

      // if requested size is too big, return an empty vector
      if( HELICES.GetSize() > s_max_sses || STRANDS.GetSize() > s_max_sses)
      {
        return large_templates;
      }

      // create an iterator to the vector of templates with the correct number of helices and strands
      storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      >::const_iterator map_itr
      (
        GetTemplates().Find( storage::Pair< size_t, size_t>( HELICES.GetSize(), STRANDS.GetSize()))
      );

      // create an iterator to the vector of templates with the correct maximum of helices and strands
      storage::Map
      <
        storage::Pair< size_t, size_t>, util::ShPtrVector< FoldTemplate>
      >::const_iterator map_itr_end
      (
        GetTemplates().Find( storage::Pair< size_t, size_t>( s_max_sses, s_max_sses))
      );

      // if the iterator is invalid
      if( map_itr == GetTemplates().End())
      {
        // set it to the beginning
        map_itr = GetTemplates().Begin();
      }
      ++map_itr_end;

      // combine the SSEs into one list
      util::SiPtrList< const SSE> combined_sses( HELICES);
      combined_sses.Append( STRANDS);

      // iterate through the map
      for( ; map_itr != map_itr_end; ++map_itr)
      {
        // iterate through these templates
        for
        (
          util::ShPtrVector< FoldTemplate>::const_iterator
            template_itr( map_itr->second.Begin()), template_itr_end( map_itr->second.End());
          template_itr != template_itr_end; ++template_itr
        )
        {
          // add the template to the vector if it has enough helices and strands
          if
          (
            map_itr->first.First() >= HELICES.GetSize() &&
            map_itr->first.Second() >= STRANDS.GetSize() &&
            map_itr->first.Second() <= s_max_sses
          )
          {
            // iterate over the SSEs
            for
            (
              util::SiPtrList< const SSE>::const_iterator sse_itr( combined_sses.Begin()),
                sse_itr_end( combined_sses.End());
              sse_itr != sse_itr_end; ++sse_itr
            )
            {
              // make sure the SSE would fit into at least one geometry
              if( !( *template_itr)->HasSimilarSizeGeometry( *sse_itr, SSE_GEOMETRY_COMPARE))
              {
                // skip to next template
                continue;
              }
            }

            // this is a suitable template so push it back
            large_templates.PushBack( *template_itr);
          }
        }
      }

      // end
      return large_templates;
    }

    //! @brief default input filename
    //! @param TEMPLATE_TYPE type of templates to use
    //! @return default input filename
    const std::string &FoldTemplateHandler::GetInputFilename( const std::string &TEMPLATE_TYPE)
    {
      // initialize with default file names
      static const std::string s_soluble_filename( score::Score::AddHistogramPath( "soluble_fold_templates.input"));
      static const std::string s_membrane_filename( score::Score::AddHistogramPath( "membrane_fold_templates.input"));
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
    storage::Set< std::string> &FoldTemplateHandler::GetExcludedPdbs()
    {
      static storage::Set< std::string> s_excluded_pdbs;

      // if a file that contains PDBs to exclude as templates was passed
      if( s_excluded_pdbs.IsEmpty() && GetFlagFoldTemplates()->GetParameterList()( 4)->GetWasSetInCommandLine())
      {
        // open the file
        io::IFStream read;
        io::File::MustOpenIFStream( read, GetFlagFoldTemplates()->GetParameterList()( 4)->GetValue());

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

  } // namespace assemble
} // namespace bcl
