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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "app/bcl_app_groups.h"
#include "app/bcl_app_interface.h"
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_pick_atom_random.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "math/bcl_math_template_instantiations.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "storage/bcl_storage_table.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeFLDResults
    //! @brief Application for analyzing fragments of molecules which are derived from a common scaffold
    //! @brief using the FocusedLibraryDesign application
    //!
    //! @author geanesar
    //! @date 05/12/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class AnalyzeFLDResults :
      public Interface
    {

    public:

    //////////
    // data //
    //////////
      util::ShPtr< command::FlagInterface> m_ScaffoldFilenameFlag;

      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      util::ShPtr< command::FlagInterface> m_PropertiesObjectFileFlag;

      util::ShPtr< command::FlagInterface> m_FragmentOutputFlag;

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      AnalyzeFLDResults();

      AnalyzeFLDResults *Clone() const
      {
        return new AnalyzeFLDResults( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_ScaffoldFilenameFlag);
        sp_cmd->AddFlag( m_InputFilenamesFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_PropertiesObjectFileFlag);
        sp_cmd->AddFlag( m_FragmentOutputFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      int Main() const
      {

        // Read in the scaffold
        io::IFStream scaffold_in;
        io::File::MustOpenIFStream( scaffold_in, m_ScaffoldFilenameFlag->GetFirstParameter()->GetValue());

        chemistry::FragmentComplete scaffold( sdf::FragmentFactory::MakeFragment( scaffold_in, sdf::e_Maintain));

        io::File::CloseClearFStream( scaffold_in);

        // Read in the list of derived molecules
        const storage::Vector< std::string> filenames( m_InputFilenamesFlag->GetStringList());
        chemistry::FragmentEnsemble derivatives;

        io::IFStream derivs_in;

        for
        (
          storage::Vector< std::string>::const_iterator itr( filenames.Begin()), itr_end( filenames.End());
          itr != itr_end;
          ++itr
        )
        {
          io::File::MustOpenIFStream( derivs_in, *itr);
          derivatives.ReadMoreFromMdl( derivs_in, sdf::e_Maintain);
          io::File::CloseClearFStream( derivs_in);
        }

        // Read in the properties object
        io::IFStream obj_input;
        io::File::MustOpenIFStream( obj_input, m_PropertiesObjectFileFlag->GetFirstParameter()->GetValue());
        const util::ObjectDataLabel properties( obj_input);
        io::File::CloseClearFStream( obj_input);

        bool write_fragments( false);
        io::OFStream output_fragment;
        if( m_FragmentOutputFlag->GetFlag() && m_FragmentOutputFlag->GetParameterList().GetSize() > 0)
        {
          write_fragments = true;
          io::File::MustOpenOFStream( output_fragment, m_FragmentOutputFlag->GetFirstParameter()->GetValue());
        }

        // Add each property to a vector
        storage::Vector< descriptor::CheminfoProperty> atom_properties;
        storage::Vector< std::string> table_header_strings( 0); // String vector used for table header

        // Add a column for molecule number, grow point, fragment number, and activity
        table_header_strings.PushBack( "Molecule_Number");
        table_header_strings.PushBack( "Grow_Point");
        table_header_strings.PushBack( "Activity");

        for
        (
          util::ObjectDataLabel::const_iterator itr( properties.Begin()), itr_end( properties.End());
          itr != itr_end;
          ++itr
        )
        {
          atom_properties.PushBack( descriptor::CheminfoProperty( *itr));
          table_header_strings.PushBack( atom_properties.Last()->GetString());
          BCL_MessageStd( "Added property: " + atom_properties.Last()->GetString());
        }

        // Set up the table used to store fragment properties
        storage::TableHeader table_header( table_header_strings);
        storage::Table< float> summed_properties_table( table_header);

        // Open the output file
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());

        // Use the CollectorValence class to collect atoms which don't have filled valences
        // in order to determine where fragments were added (grow points)
        const util::SiPtrList< const chemistry::AtomConformationalInterface>
        grow_points_list( chemistry::CollectorValence().Collect( scaffold));

        // Convert the grow point list to a vector for easier handling later
        storage::Vector< size_t> grow_points_vector;
        for
        (
          util::SiPtrList< const chemistry::AtomConformationalInterface>::const_iterator
            itr( grow_points_list.Begin()), itr_end( grow_points_list.End());
          itr != itr_end;
          ++itr
        )
        {
          grow_points_vector.PushBack( scaffold.GetAtomIndex( **itr));
        }

        // Graph maker, self explanatory
        chemistry::ConformationGraphConverter conformation_graph_maker
        (
          chemistry::ConformationGraphConverter::e_ElementType,
          chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
        );

        // create a graph for the scaffold
        graph::ConstGraph< size_t, size_t> scaffold_graph( conformation_graph_maker( scaffold));

        // Iterate through each molecule (and calculate properties on each fragment of the molecule)
        size_t molecule_index( 0);
        for
        (
          chemistry::FragmentEnsemble::iterator itr( derivatives.Begin()), itr_end( derivatives.End());
          itr != itr_end;
          ++itr, ++molecule_index
        )
        {

          // Vector to store atom-wise molecule properties in
          storage::Vector< linal::Vector< float> > property_vector;

          // Run SetObject on each property using our moleculebefore doing anything else,
          // then calculate the properties for the molecule. We will need this later for summing atom properties
          for
          (
            storage::Vector< descriptor::CheminfoProperty>::iterator itr_prop( atom_properties.Begin()),
              itr_prop_end( atom_properties.End());
            itr_prop != itr_prop_end;
            ++itr_prop
          )
          {
            ( *itr_prop)->SetObject( *itr);
            property_vector.PushBack( ( *itr_prop)->CollectValuesOnEachElementOfObject( *itr));
          }

          // Create a graph of this molecule
          graph::ConstGraph< size_t, size_t> mol_graph( conformation_graph_maker( *itr));

          graph::ConstGraph< util::SiPtr< const chemistry::AtomConformationalInterface>, size_t> mol_atom_graph
          (
            conformation_graph_maker.CreateGraphWithAtoms( *itr)
          );

          // Find the isomorphism between this molecule and the scaffold in order to find grow points
          graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
          isomorphism.SetGraphExternalOwnership( mol_graph);
          isomorphism.SetSubgraphExternalOwnership( scaffold_graph);

          // If an isomorphism doesn't exist, this scaffold doesn't exist in the molecule, skip it
          if( !isomorphism.FindIsomorphism())
          {
            BCL_MessageStd( "Scaffold is not present in derivative.  Skipping.");
            continue;
          }

          // We need the indices of the *molecule* that correspond to the grow points on the scaffold
          // Use the isomorphism to find out what these indices are
          const storage::Vector< size_t> scaffold_to_mol_iso( isomorphism.GetIsomorphism());

          // Fill in the vector of grow point indices using the calculated isomorphism
          // Iterate through the original grow points vector (grow points for the scaffold) and
          // use this as the index into the isomorphism
          // This vector holds indices from the current molecule which correspond to grow points
          storage::Vector< size_t> grow_point_indices;
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_growpt( grow_points_vector.Begin()),
              itr_growpt_end( grow_points_vector.End());
            itr_growpt != itr_growpt_end;
            ++itr_growpt
          )
          {
            grow_point_indices.PushBack( scaffold_to_mol_iso( *itr_growpt));
          }

          BCL_MessageVrb( "Number of grow points: " + util::Format()( grow_points_vector.GetSize()));

          BCL_MessageDbg( "grow point indices: " + util::Format()( grow_point_indices));

          // 1) Find edges that connect the scaffold to the larger molecule;
          // 2) Remove the found edges, then fetch the connected pieces associated with the connected atom
          // 3) Store the connected pieces to calculate properties later
          storage::Vector< size_t> neighboring_atoms; // Neighboring atoms to the current grow point

          storage::Vector< storage::Pair< size_t, storage::Vector< size_t> > > fragments; // This stores fragments (as a list of indices), and grow points each is associated with

          // Iterate over grow points:
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_grow_point( grow_point_indices.Begin()), itr_grow_point_end( grow_point_indices.End());
            itr_grow_point != itr_grow_point_end;
            ++itr_grow_point
          )
          {
            BCL_MessageVrb( "grow_index: " + util::Format()( *itr_grow_point));

            // Now:
            // Get indices of neighboring vertices, then check if they belong to the scaffold.  If they do not,
            // remove the edge between the grow point and the neighbor.  Store the neighbor for later queries

            // Holds the neighboring indices of this grow point
            storage::Vector< size_t> neighbor_indices( mol_graph.GetNeighborIndices( *itr_grow_point));

            // Holds neighbor indices associated with fragments (i.e. non-scaffold parts of the molecule)
            storage::Vector< size_t> fragment_indices;

            for
            (
              storage::Vector< size_t>::const_iterator itr_neighbor( neighbor_indices.Begin()),
                itr_neighbor_end( neighbor_indices.End());
              itr_neighbor != itr_neighbor_end;
              ++itr_neighbor
            )
            {
              // If the grow point (*itr) is not in the isomorphism (scaffold_to_mol_iso), then it is not part of the
              // scaffold.  Therefore, remove the edge and keep track of the neighbor index
              if( scaffold_to_mol_iso.Find( *itr_neighbor) >= scaffold_to_mol_iso.GetSize())
              {
                BCL_MessageVrb( "neighbor index: " + util::Format()( *itr) + " to grow index: "
                  + util::Format()( *itr_grow_point));

                fragment_indices.PushBack( *itr_neighbor);

                mol_graph.RemoveEdge( *itr_grow_point, *itr_neighbor);
              }
            }

            // Now, all edges to the fragments attached at this grow point are removed, so we can collect the
            // connected pieces of the graph associated with each (formerly attached) neighbor
            for
            (
              storage::Vector< size_t>::iterator itr_frag_index( fragment_indices.Begin()), itr_frag_index_end( fragment_indices.End());
              itr_frag_index != itr_frag_index_end;
              ++itr_frag_index
            )
            {
              // The the fragment indices
              util::ShPtr< storage::Vector< size_t> > sp_this_fragment_indices
              (
                graph::Connectivity::GetVerticesReachableFrom
                (
                  mol_graph,
                  *itr_frag_index
                )
              );

              // Add the fragment to the fragment list using the grow point from the scaffold and the list of indices
              // that we found above
              fragments.PushBack
              (
                storage::Pair< size_t, storage::Vector< size_t> >
                (
                  scaffold_to_mol_iso.Find( *itr_grow_point),
                  *sp_this_fragment_indices
                )
              );
            }
          } // End of iterating over grow points

          // The graph has been disconnected to isolate fragments from the scaffold.  Now we want to get the parts of
          // the graph that are still connected (stored as vectors of indices)
          storage::List< storage::Vector< size_t> > connected_parts( graph::Connectivity::GetComponents( mol_graph));

          // Iterate through each fragment and add up properties for the atoms
          for
          (
            storage::Vector< storage::Pair< size_t, storage::Vector< size_t> > >::iterator
              itr_fragment( fragments.Begin()), itr_fragment_end( fragments.End());
            itr_fragment != itr_fragment_end;
            ++itr_fragment
          )
          {

            // If we want to write out the fragments that were found, do so
            if( write_fragments)
            {
              chemistry::FragmentComplete fragment
              (
                conformation_graph_maker.CreateAtomsFromGraph( mol_atom_graph.GetSubgraph( itr_fragment->Second())),
                ""
              );
              fragment.WriteMDL( output_fragment);
            }

            // Vector which will hold the summed properties for this fragment
            storage::Vector< float> this_fragment_summed_properties( 0);

            // Row name we will store in the table
            std::string row_name = "mol" + util::Format()( molecule_index) + "_" + itr->GetSumFormula() +
                "_gp" + util::Format()( itr_fragment->First());

            // Add molecule number, grow point, fragment number, and activity to the table
            this_fragment_summed_properties.PushBack( molecule_index);
            this_fragment_summed_properties.PushBack( itr_fragment->First());

            std::string mol_activity( itr->GetMDLProperty( "Activity"));
            if( mol_activity.length() == 0)
            {
              this_fragment_summed_properties.PushBack( 0.0);
            }
            else
            {
              this_fragment_summed_properties.PushBack
              (
                util::ConvertStringToNumericalValue< float>( mol_activity)
              );
            }

            // Iterate through each property to collect the values for each atom in this fragment
            for
            (
              storage::Vector< linal::Vector< float> >::iterator itr_prop( property_vector.Begin()),
                itr_prop_end( property_vector.End());
              itr_prop != itr_prop_end;
              ++itr_prop
            )
            {
              float result( 0.0);

              // Iterate across each atom of the fragment and add the property value for that atom to result
              for
              (
                storage::Vector< size_t>::iterator itr_atom_index( itr_fragment->Second().Begin()),
                  itr_atom_index_end( itr_fragment->Second().End());
                itr_atom_index != itr_atom_index_end;
                ++itr_atom_index
              )
              {
                // Double check to make sure we are accessing a properly sized vector
                if( ( *itr_atom_index) >= itr_prop->GetSize())
                {
                  BCL_MessageStd( util::Format()( *itr_atom_index) + " is bigger than " + util::Format()( itr_prop->GetSize()));
                  BCL_Exit( "Property is not an atom-wise property! Exitting", -1);
                }
                result += itr_prop->operator()( *itr_atom_index);
              }
              this_fragment_summed_properties.PushBack( result);
            } // End of iterating through properties

            // We have now summed up the properties for this fragment.  Now, we need to check to see if there is
            // a row with the name we want already present.  If there is, get the data from it, add the atom properties
            // to the ones we have calculated, and then remove the original from the table
            if( summed_properties_table.HasRow( row_name))
            {
              storage::Vector< float> original_props( summed_properties_table[ row_name].GetData());

              // Get the column which contains "Activity", since this will be the last row that won't change based
              // on the grow point (molecule-general?)
              size_t activity_column( summed_properties_table.GetHeader()[ "Activity"]);

              // Get an iterator for the original property vector
              storage::Vector< float>::iterator itr_orig_prop( original_props[ activity_column + 1]);
              for
              (
                storage::Vector< float>::iterator itr_new_prop( this_fragment_summed_properties[ activity_column + 1]),
                  itr_new_prop_end( this_fragment_summed_properties.End());
                itr_new_prop != itr_new_prop_end;
                ++itr_new_prop, ++itr_orig_prop
              )
              {
                // Add the old value to the new value
                ( *itr_new_prop) += ( *itr_orig_prop);
              }

              // Remove the original row from the table so that the new one can be inserted
              summed_properties_table.RemoveRow( row_name);
            }

            // Insert the new row into the table
            summed_properties_table.InsertRow
            (
              row_name,
              this_fragment_summed_properties
            );
          }
        }

        if( write_fragments)
        {
          io::File::CloseClearFStream( output_fragment);
        }

        summed_properties_table.WriteFormatted( output);
        BCL_MessageStd( "Wrote table to output");
        io::File::CloseClearFStream( output);
        return 0;
      } // Main()

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

    private:

      static const ApplicationType AnalyzeFLDResults_Instance;

    }; // AnalyzeFLDResults

      //! @brief standard constructor
    AnalyzeFLDResults::AnalyzeFLDResults() :
      m_ScaffoldFilenameFlag
      (
        new command::FlagStatic
        (
          "scaffold", "filename for the scaffold for the molecules",
          command::Parameter
          (
            "filename",
            "filename for scaffold sdf",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename", "flag selecting the output file name",
          command::Parameter
          (
            "filename", "filename for output sdf"
          )
        )
      ),
      m_InputFilenamesFlag
      (
        new command::FlagDynamic
        (
          "input_filenames",
          "filenames for input sdf of derived structures",
          command::Parameter
          (
            "filenames of derived structure sdfs",
            "name of files containing derived structures",
            command::ParameterCheckFileExistence()
          ),
          1,
          21
        )
      ),
      m_PropertiesObjectFileFlag
      (
        new command::FlagStatic
        (
          "code_object", "filename containing code object for atom properties",
          command::Parameter
          (
            "filename",
            "filename for code object",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_FragmentOutputFlag
      (
        new command::FlagDynamic
        (
          "output_fragments",
          "filename to write out fragments",
          command::Parameter
          (
            "output_fragments",
            "name of the file the fragments will be written to"
          ),
          0,
          1
        )
      )
    {
    }

    const ApplicationType AnalyzeFLDResults::AnalyzeFLDResults_Instance
    (
      GetAppGroups().AddAppToGroup( new AnalyzeFLDResults(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
