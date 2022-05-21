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
#include "nmr/bcl_nmr_star_tags.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StarTags::StarTags() :
      e_DistTreeConstraintID(              AddEnum( "DistTreeConstraintID",              StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Constraint_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistTreeNodeID(                    AddEnum( "DistTreeNodeID",                    StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Node_ID",                        util::CPPDataTypes::e_SizeT))),
      e_DistTreeDownNodeID(                AddEnum( "DistTreeDownNodeID",                StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Down_node_ID",                   util::CPPDataTypes::e_SizeT))),
      e_DistTreeRightNodeID(               AddEnum( "DistTreeRightNodeID",               StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Right_node_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistTreeLogicOperation(            AddEnum( "DistTreeLogicOperation",            StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Logic_operation",                util::CPPDataTypes::e_String))),
      e_DistTreeEntryID(                   AddEnum( "DistTreeEntryID",                   StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistTreeDistanceConstraintListID(  AddEnum( "DistTreeDistanceConstraintListID",  StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Distance_constraint_list_ID",    util::CPPDataTypes::e_SizeT))),

      e_DistTreeNodeMemberConstraintID(    AddEnum( "DistTreeNodeMemberConstraintID",    StarTagData( GetStarTagCategories().e_DistConstraint,      "Tree_node_member_constraint_ID", util::CPPDataTypes::e_SizeT))),
      e_DistTreeNodeMemberNodeID(          AddEnum( "DistTreeNodeMemberNodeID",          StarTagData( GetStarTagCategories().e_DistConstraint,      "Tree_node_member_node_ID",       util::CPPDataTypes::e_SizeT))),
      e_DistConstraintTreeNodeMemberID(    AddEnum( "DistConstraintTreeNodeMemberID",    StarTagData( GetStarTagCategories().e_DistConstraint,      "Constraint_tree_node_member_ID", util::CPPDataTypes::e_SizeT))),
      e_DistAssemblyAtomID(                AddEnum( "DistAssemblyAtomID",                StarTagData( GetStarTagCategories().e_DistConstraint,      "Assemble_atom_ID",               util::CPPDataTypes::e_SizeT))),
      e_DistEntityAssemblyID(              AddEnum( "DistEntityAssemblyID",              StarTagData( GetStarTagCategories().e_DistConstraint,      "Entity_assembly_ID",             util::CPPDataTypes::e_SizeT))),
      e_DistEntityID(                      AddEnum( "DistEntityID",                      StarTagData( GetStarTagCategories().e_DistConstraint,      "Entity_ID",                      util::CPPDataTypes::e_SizeT))),
      e_DistCompIndexID(                   AddEnum( "DistCompIndexID",                   StarTagData( GetStarTagCategories().e_DistConstraint,      "Comp_index_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistSeqID(                         AddEnum( "DistSeqID",                         StarTagData( GetStarTagCategories().e_DistConstraint,      "Seq_ID",                         util::CPPDataTypes::e_SizeT))),
      e_DistCompID(                        AddEnum( "DistCompID",                        StarTagData( GetStarTagCategories().e_DistConstraint,      "Comp_ID",                        util::CPPDataTypes::e_String))),
      e_DistAtomID(                        AddEnum( "DistAtomID",                        StarTagData( GetStarTagCategories().e_DistConstraint,      "Atom_ID",                        util::CPPDataTypes::e_String))),
      e_DistResonanceID(                   AddEnum( "DistResonanceID",                   StarTagData( GetStarTagCategories().e_DistConstraint,      "Resonance_ID",                   util::CPPDataTypes::e_SizeT))),
      e_DistAuthAsymID(                    AddEnum( "DistAuthAsymID",                    StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_asym_ID",                   util::CPPDataTypes::e_String))),
      e_DistAuthSeqID(                     AddEnum( "DistAuthSeqID",                     StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_seq_ID",                    util::CPPDataTypes::e_String))),
      e_DistAuthCompID(                    AddEnum( "DistAuthCompID",                    StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_comp_ID",                   util::CPPDataTypes::e_String))),
      e_DistAuthAtomID(                    AddEnum( "DistAuthAtomID",                    StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_atom_ID",                   util::CPPDataTypes::e_String))),
      e_DistEntryID(                       AddEnum( "DistEntryID",                       StarTagData( GetStarTagCategories().e_DistConstraint,      "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistDistanceConstraintListID(      AddEnum( "DistDistanceConstraintListID",      StarTagData( GetStarTagCategories().e_DistConstraint,      "Distance_constraint_list_ID",    util::CPPDataTypes::e_SizeT))),

      e_DistListSfCatetory(                AddEnum( "DistListSfCatetory",                StarTagData( GetStarTagCategories().e_DistConstraintList,  "Sf_category",                    util::CPPDataTypes::e_String))),
      e_DistListEntryID(                   AddEnum( "DistListEntryID",                   StarTagData( GetStarTagCategories().e_DistConstraintList,  "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistListSfID(                      AddEnum( "DistListSfID",                      StarTagData( GetStarTagCategories().e_DistConstraintList,  "Sf_ID",                          util::CPPDataTypes::e_SizeT))),
      e_DistListID(                        AddEnum( "DistListID",                        StarTagData( GetStarTagCategories().e_DistConstraintList,  "ID",                             util::CPPDataTypes::e_SizeT))),
      e_DistListConstraintType(            AddEnum( "DistListConstraintType",            StarTagData( GetStarTagCategories().e_DistConstraintList,  "Constraint_type",                util::CPPDataTypes::e_String))),
      e_DistListConstraintFileID(          AddEnum( "DistListConstraintFileID",          StarTagData( GetStarTagCategories().e_DistConstraintList,  "Constraint_file_ID",             util::CPPDataTypes::e_SizeT))),
      e_DistListBlockID(                   AddEnum( "DistListBlockID",                   StarTagData( GetStarTagCategories().e_DistConstraintList,  "Block_ID",                       util::CPPDataTypes::e_SizeT))),
      e_DistListDetails(                   AddEnum( "DistListDetails",                   StarTagData( GetStarTagCategories().e_DistConstraintList,  "Details",                        util::CPPDataTypes::e_String))),

      e_DistValueConstraintID(             AddEnum( "DistValueConstraintID",             StarTagData( GetStarTagCategories().e_DistConstraintValue, "Constraint_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistValueTreeNodeID(               AddEnum( "DistValueTreeNodeID",               StarTagData( GetStarTagCategories().e_DistConstraintValue, "Tree_node_ID",                   util::CPPDataTypes::e_SizeT))),
      e_DistValueSourceExperimentID(       AddEnum( "DistValueSourceExperimentID",       StarTagData( GetStarTagCategories().e_DistConstraintValue, "Source_experiment_ID",           util::CPPDataTypes::e_SizeT))),
      e_DistValueSpectralPeakID(           AddEnum( "DistValueSpectralPeakID",           StarTagData( GetStarTagCategories().e_DistConstraintValue, "Spectral_peak_ID",               util::CPPDataTypes::e_SizeT))),
      e_DistValueIntensityVal(             AddEnum( "DistValueIntensityVal",             StarTagData( GetStarTagCategories().e_DistConstraintValue, "Intensity_val",                  util::CPPDataTypes::e_Double))),
      e_DistValueIntensityLowerValErr(     AddEnum( "DistValueIntensityLowerValErr",     StarTagData( GetStarTagCategories().e_DistConstraintValue, "Intensity_lower_val_err",        util::CPPDataTypes::e_Double))),
      e_DistValueIntensityUpperValErr(     AddEnum( "DistValueIntensityUpperValErr",     StarTagData( GetStarTagCategories().e_DistConstraintValue, "Intensity_upper_val_err",        util::CPPDataTypes::e_Double))),
      e_DistValueDistanceVal(              AddEnum( "DistValueDistanceVal",              StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_val",                   util::CPPDataTypes::e_Double))),
      e_DistValueDistanceLowerBoundVal(    AddEnum( "DistValueDistanceLowerBoundVal",    StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_lower_bound_val",       util::CPPDataTypes::e_Double))),
      e_DistValueDistanceUpperBoundVal(    AddEnum( "DistValueDistanceUpperBoundVal",    StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_upper_bound_val",       util::CPPDataTypes::e_Double))),
      e_DistValueEntryID(                  AddEnum( "DistValueEntryID",                  StarTagData( GetStarTagCategories().e_DistConstraintValue, "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistValueDistanceConstraintListID( AddEnum( "DistValueDistanceConstraintListID", StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_constraint_list_ID",    util::CPPDataTypes::e_SizeT))),

      e_RDCID(                             AddEnum( "RDCID",                             StarTagData( GetStarTagCategories().e_RDCConstraint,       "ID",                             util::CPPDataTypes::e_SizeT))),
      e_RDCAssemblyAtomID1(                AddEnum( "RDCAssemblyAtomID1",                StarTagData( GetStarTagCategories().e_RDCConstraint,       "Assembly_atom_ID_1",             util::CPPDataTypes::e_SizeT))),
      e_RDCEntityAssemblyID1(              AddEnum( "RDCEntityAssemblyID1",              StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_assembly_ID_1",           util::CPPDataTypes::e_SizeT))),
      e_RDCEntityID1(                      AddEnum( "RDCEntityID1",                      StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_ID_1",                    util::CPPDataTypes::e_SizeT))),
      e_RDCCompIndexID1(                   AddEnum( "RDCCompIndexID1",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_index_ID_1",                util::CPPDataTypes::e_SizeT))),
      e_RDCSeqID1(                         AddEnum( "RDCSeqID1",                         StarTagData( GetStarTagCategories().e_RDCConstraint,       "Seq_ID_1",                       util::CPPDataTypes::e_SizeT))),
      e_RDCCompID1(                        AddEnum( "RDCCompID1",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_ID_1",                      util::CPPDataTypes::e_String))),
      e_RDCAtomID1(                        AddEnum( "RDCAtomID1",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Atom_ID_1",                      util::CPPDataTypes::e_String))),
      e_RDCResonanceID1(                   AddEnum( "RDCResonanceID1",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Resonance_ID_1",                 util::CPPDataTypes::e_SizeT))),
      e_RDCAssemblyAtomID2(                AddEnum( "RDCAssemblyAtomID2",                StarTagData( GetStarTagCategories().e_RDCConstraint,       "Assembly_atom_ID_2",             util::CPPDataTypes::e_SizeT))),
      e_RDCEntityAssemblyID2(              AddEnum( "RDCEntityAssemblyID2",              StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_assembly_ID_2",           util::CPPDataTypes::e_SizeT))),
      e_RDCEntityID2(                      AddEnum( "RDCEntityID2",                      StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_ID_2",                    util::CPPDataTypes::e_SizeT))),
      e_RDCCompIndexID2(                   AddEnum( "RDCCompIndexID2",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_index_ID_2",                util::CPPDataTypes::e_SizeT))),
      e_RDCSeqID2(                         AddEnum( "RDCSeqID2",                         StarTagData( GetStarTagCategories().e_RDCConstraint,       "Seq_ID_2",                       util::CPPDataTypes::e_SizeT))),
      e_RDCCompID2(                        AddEnum( "RDCCompID2",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_ID_2",                      util::CPPDataTypes::e_String))),
      e_RDCAtomID2(                        AddEnum( "RDCAtomID2",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Atom_ID_2",                      util::CPPDataTypes::e_String))),
      e_RDCResonanceID2(                   AddEnum( "RDCResonanceID2",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Resonance_ID_2",                 util::CPPDataTypes::e_SizeT))),
      e_RDCVal(                            AddEnum( "RDCVal",                            StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_val",                        util::CPPDataTypes::e_Double))),
      e_RDCLowerBound(                     AddEnum( "RDCLowerBound",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_lower_bound",                util::CPPDataTypes::e_Double))),
      e_RDCUpperBound(                     AddEnum( "RDCUpperBound",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_upper_bound",                util::CPPDataTypes::e_Double))),
      e_RDCValErr(                         AddEnum( "RDCValErr",                         StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_val_err",                    util::CPPDataTypes::e_Double))),
      e_RDCBondLength(                     AddEnum( "RDCBondLength",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_bond_length",                util::CPPDataTypes::e_Double))),
      e_RDCSourceExperimentID(             AddEnum( "RDCSourceExperimentID",             StarTagData( GetStarTagCategories().e_RDCConstraint,       "Source_experiment_ID",           util::CPPDataTypes::e_SizeT))),
      e_RDCAuthAsymID1(                    AddEnum( "RDCAuthAsymID1",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_asym_ID_1",                 util::CPPDataTypes::e_String))),
      e_RDCAuthSeqID1(                     AddEnum( "RDCAuthSeqID1",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_seq_ID_1",                  util::CPPDataTypes::e_String))),
      e_RDCAuthCompID1(                    AddEnum( "RDCAuthCompID1",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_comp_ID_1",                 util::CPPDataTypes::e_String))),
      e_RDCAuthAtomID1(                    AddEnum( "RDCAuthAtomID1",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_atom_ID_1",                 util::CPPDataTypes::e_String))),
      e_RDCAuthAsymID2(                    AddEnum( "RDCAuthAsymID2",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_asym_ID_2",                 util::CPPDataTypes::e_String))),
      e_RDCAuthSeqID2(                     AddEnum( "RDCAuthSeqID2",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_seq_ID_2",                  util::CPPDataTypes::e_String))),
      e_RDCAuthCompID2(                    AddEnum( "RDCAuthCompID2",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_comp_ID_2",                 util::CPPDataTypes::e_String))),
      e_RDCAuthAtomID2(                    AddEnum( "RDCAuthAtomID2",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_atom_ID_2",                 util::CPPDataTypes::e_String))),
      e_RDCEntryID(                        AddEnum( "RDCEntryID",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_RDCRDCConstraintListID(            AddEnum( "RDCRDCConstraintListID",            StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_constraint_list_ID",         util::CPPDataTypes::e_SizeT))),

      e_RDCListSfCategory(                 AddEnum( "RDCListSfCategory",                 StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Sf_category",                    util::CPPDataTypes::e_String))),
      e_RDCListEntryID(                    AddEnum( "RDCListEntryID",                    StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_RDCListID(                         AddEnum( "RDCListID",                         StarTagData( GetStarTagCategories().e_RDCConstraintList,   "ID",                             util::CPPDataTypes::e_SizeT))),
      e_RDCListConstraintFileID(           AddEnum( "RDCListConstraintFileID",           StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Constraint_file_ID",             util::CPPDataTypes::e_SizeT))),
      e_RDCListBlockID(                    AddEnum( "RDCListBlockID",                    StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Block_ID",                       util::CPPDataTypes::e_SizeT))),
      e_RDCListDetais(                     AddEnum( "RDCListDetais",                     StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Details",                        util::CPPDataTypes::e_String)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTags::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a tag from a string
    //! @param CATEGORY this tag belongs to
    //! @param DESCRIPTION string description for the enum
    //! @return a tag from a string
    StarTag StarTags::GetTagFromString( const StarTagCategory &CATEGORY, const std::string &DESCRIPTION)
    {
      // iterate through the tag categories
      for
      (
        StarTags::const_iterator tag_itr( GetStarTags().Begin()), tag_itr_end( GetStarTags().End());
        tag_itr != tag_itr_end;
        ++tag_itr
      )
      {
        // if the category has the same category and description
        if( ( *tag_itr)->GetTagCategory() == CATEGORY && ( *tag_itr)->GetDescription() == DESCRIPTION)
        {
          // return it
          return *tag_itr;
        }
      }

      // no matching description found
      return GetStarTags().e_Undefined;
    }

    //! @brief read Star file from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    storage::Map
    <
      StarTagCategory,
      storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
    > StarTags::ReadStarFile( std::istream &ISTREAM)
    {
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > line_data;

      // read through the file
      while( !ISTREAM.eof())
      {
        // store the line
        std::string buffer;
        std::getline( ISTREAM, buffer);

        // check if it is a start of a data loop
        if( util::TrimString( buffer) == "loop_")
        {
          // get the next line
          std::getline( ISTREAM, buffer);

          // check if the line contains a tag category-tag pair of the format "_Category.Tag"
          storage::Vector< std::string> split_line( util::SplitString( buffer, "."));
          StarTagCategory category
          (
            StarTagCategories::GetCategoryFromString( util::TrimString( split_line.FirstElement()))
          );

          // true if the category is defined
          if( category != GetStarTagCategories().e_Undefined)
          {
            // read in the tag associated with the valid category
            StarTag tag( GetTagFromString( category, util::TrimString( split_line.LastElement())));

            // add the first tag to the vector in the map
            line_data[ category].First().PushBack( tag);

            // read the rest of the tags
            while( !ISTREAM.eof())
            {
              // read in the next line
              std::getline( ISTREAM, buffer);
              split_line = util::SplitString( buffer, ".");

              // if the split line does not have 2 elements, there are no more tags to read in so break
              if( split_line.GetSize() != 2)
              {
                break;
              }

              // add the tag
              line_data[ category].First().PushBack
              (
                GetTagFromString( category, util::TrimString( split_line.LastElement()))
              );
            }

            // iterate through any blank lines
            while( !ISTREAM.eof() && util::TrimString( buffer).empty())
            {
              std::getline( ISTREAM, buffer);
            }

            // iterate through the rows of data that correspond to the tags just read in
            while( !ISTREAM.eof() && util::TrimString( buffer) != "stop_")
            {
              // add the line to the list of strings
              line_data[ category].Second().PushBack( util::TrimString( buffer));
              std::getline( ISTREAM, buffer);
            }
          }

        } // end of data loop
      } // eof
      // end
      return line_data;
    }

    //! @brief gets a chain id from a string since the star format allows authors to identify chains however they want
    //! @param CURRENT_ID current string to be converted
    //! @param PREVIOUS_IDS previously seen strings mapped to chain ids
    //! @return chain id
    char StarTags::GetChainID
    (
      const std::string &CURRENT_ID,
      storage::Map< std::string, char> &PREVIOUS_IDS
    )
    {
      // check the map for the id
      const storage::Map< std::string, char>::const_iterator find_itr( PREVIOUS_IDS.Find( CURRENT_ID));

      // if it is found
      if( find_itr != PREVIOUS_IDS.End())
      {
        // return the corresponding chain
        return find_itr->second;
      }

      // it was not found, so check if the string is 1 character
      if( CURRENT_ID.size() == 1 && CURRENT_ID[ 0] >= 'A' && CURRENT_ID[ 0] <= 'Z')
      {
        // add to map and return
        const char chain_id( CURRENT_ID[ 0]);
        PREVIOUS_IDS[ CURRENT_ID] = chain_id;
        return chain_id;
      }

      // a new, 2+ character id has been passed
      // check if any ids have been set
      if( PREVIOUS_IDS.IsEmpty())
      {
        // set this id to A
        const char initial_id( 'A');
        PREVIOUS_IDS[ CURRENT_ID] = initial_id;
        return initial_id;
      }

      // create a set to find the largets current chain id char
      storage::Set< char> previous_chain_ids;

      // iterate over the map
      for
      (
        storage::Map< std::string, char>::const_iterator map_itr( PREVIOUS_IDS.Begin()),
          map_itr_end( PREVIOUS_IDS.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        previous_chain_ids.Insert( map_itr->second);
      }

      // set the id and return
      char max_chain_id( *previous_chain_ids.ReverseBegin());
      ++max_chain_id;
      PREVIOUS_IDS[ CURRENT_ID] = max_chain_id;
      return max_chain_id;
    }

    //! @brief get access to all star tags
    const StarTags &GetStarTags()
    {
      return StarTags::GetEnums();
    }

  } // namespace nmr

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< nmr::StarTagData, nmr::StarTags>;

  } // namespace util
} // namespace bcl
