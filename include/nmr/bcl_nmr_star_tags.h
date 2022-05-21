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

#ifndef BCL_NMR_STAR_TAGS_H_
#define BCL_NMR_STAR_TAGS_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_nmr_star_tag_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StarTags
    //! @brief Enum-based class that stores NMR-Star tags to be used when reading in STAR files using the
    //!        handler.  See http://www.bmrb.wisc.edu/dictionary/3.1html/SuperGroupPage.html for more info.
    //!
    //! @see @link example_nmr_star_tags.cpp @endlink
    //! @author weinerbe
    //! @date Jun 22, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StarTags :
      public util::Enumerate< StarTagData, StarTags>
    {
      friend class util::Enumerate< StarTagData, StarTags>;

    private:

    //////////
    // data //
    //////////

    public:

    //////////
    // data //
    //////////

      //! distance tree tags
      const StarTag e_DistTreeConstraintID;              //!< e_DistTreeConstraintID;
      const StarTag e_DistTreeNodeID;                    //!< e_DistTreeNodeID;
      const StarTag e_DistTreeDownNodeID;                //!< e_DistTreeDownNodeID;
      const StarTag e_DistTreeRightNodeID;               //!< e_DistTreeRightNodeID;
      const StarTag e_DistTreeLogicOperation;            //!< e_DistTreeLogicOperation;
      const StarTag e_DistTreeEntryID;                   //!< e_DistTreeEntryID;
      const StarTag e_DistTreeDistanceConstraintListID;  //!< e_DistTreeDistanceConstraintListID;

      //! distance constraint tags
      const StarTag e_DistTreeNodeMemberConstraintID;    //!< e_DistTreeNodeMemberConstraintID;
      const StarTag e_DistTreeNodeMemberNodeID;          //!< e_DistTreeNodeMemberNodeID;
      const StarTag e_DistConstraintTreeNodeMemberID;    //!< e_DistConstraintTreeNodeMemberID;
      const StarTag e_DistAssemblyAtomID;                //!< e_DistAssemblyAtomID;
      const StarTag e_DistEntityAssemblyID;              //!< e_DistEntityAssemblyID;
      const StarTag e_DistEntityID;                      //!< e_DistEntityID;
      const StarTag e_DistCompIndexID;                   //!< e_DistCompIndexID;
      const StarTag e_DistSeqID;                         //!< e_DistSeqID;
      const StarTag e_DistCompID;                        //!< e_DistCompID;
      const StarTag e_DistAtomID;                        //!< e_DistAtomID;
      const StarTag e_DistResonanceID;                   //!< e_DistResonanceID;
      const StarTag e_DistAuthAsymID;                    //!< e_DistAuthAsymID;
      const StarTag e_DistAuthSeqID;                     //!< e_DistAuthSeqID;
      const StarTag e_DistAuthCompID;                    //!< e_DistAuthCompID;
      const StarTag e_DistAuthAtomID;                    //!< e_DistAuthAtomID;
      const StarTag e_DistEntryID;                       //!< e_DistEntryID;
      const StarTag e_DistDistanceConstraintListID;      //!< e_DistDistanceConstraintListID;

      //! distance list tags
      const StarTag e_DistListSfCatetory;                //!< e_DistListSfCatetory;
      const StarTag e_DistListEntryID;                   //!< e_DistListEntryID;
      const StarTag e_DistListSfID;                      //!< e_DistListSfID;
      const StarTag e_DistListID;                        //!< e_DistListID;
      const StarTag e_DistListConstraintType;            //!< e_DistListConstraintType;
      const StarTag e_DistListConstraintFileID;          //!< e_DistListConstraintFileID;
      const StarTag e_DistListBlockID;                   //!< e_DistListBlockID;
      const StarTag e_DistListDetails;                   //!< e_DistListDetails;

      //! distance value tags
      const StarTag e_DistValueConstraintID;             //!< e_DistValueConstraintID;
      const StarTag e_DistValueTreeNodeID;               //!< e_DistValueTreeNodeID;
      const StarTag e_DistValueSourceExperimentID;       //!< e_DistValueSourceExperimentID;
      const StarTag e_DistValueSpectralPeakID;           //!< e_DistValueSpectralPeakID;
      const StarTag e_DistValueIntensityVal;             //!< e_DistValueIntensityVal;
      const StarTag e_DistValueIntensityLowerValErr;     //!< e_DistValueIntensityLowerValErr;
      const StarTag e_DistValueIntensityUpperValErr;     //!< e_DistValueIntensityUpperValErr;
      const StarTag e_DistValueDistanceVal;              //!< e_DistValueDistanceVal;
      const StarTag e_DistValueDistanceLowerBoundVal;    //!< e_DistValueDistanceLowerBoundVal;
      const StarTag e_DistValueDistanceUpperBoundVal;    //!< e_DistValueDistanceUpperBoundVal;
      const StarTag e_DistValueEntryID;                  //!< e_DistValueEntryID;
      const StarTag e_DistValueDistanceConstraintListID; //!< e_DistValueDistanceConstraintListID;

      //! RDC constraint tags
      const StarTag e_RDCID;                             //!< e_RDCID;
      const StarTag e_RDCAssemblyAtomID1;                //!< e_RDCAssemblyAtomID1;
      const StarTag e_RDCEntityAssemblyID1;              //!< e_RDCEntityAssemblyID1;
      const StarTag e_RDCEntityID1;                      //!< e_RDCEntityID1;
      const StarTag e_RDCCompIndexID1;                   //!< e_RDCCompIndexID1;
      const StarTag e_RDCSeqID1;                         //!< e_RDCSeqID1;
      const StarTag e_RDCCompID1;                        //!< e_RDCCompID1;
      const StarTag e_RDCAtomID1;                        //!< e_RDCAtomID1;
      const StarTag e_RDCResonanceID1;                   //!< e_RDCResonanceID1;
      const StarTag e_RDCAssemblyAtomID2;                //!< e_RDCAssemblyAtomID2;
      const StarTag e_RDCEntityAssemblyID2;              //!< e_RDCEntityAssemblyID2;
      const StarTag e_RDCEntityID2;                      //!< e_RDCEntityID2;
      const StarTag e_RDCCompIndexID2;                   //!< e_RDCCompIndexID2;
      const StarTag e_RDCSeqID2;                         //!< e_RDCSeqID2;
      const StarTag e_RDCCompID2;                        //!< e_RDCCompID2;
      const StarTag e_RDCAtomID2;                        //!< e_RDCAtomID2;
      const StarTag e_RDCResonanceID2;                   //!< e_RDCResonanceID2;
      const StarTag e_RDCVal;                            //!< e_RDCVal;
      const StarTag e_RDCLowerBound;                     //!< e_RDCLowerBound;
      const StarTag e_RDCUpperBound;                     //!< e_RDCUpperBound;
      const StarTag e_RDCValErr;                         //!< e_RDCValErr;
      const StarTag e_RDCBondLength;                     //!< e_RDCBondLength;
      const StarTag e_RDCSourceExperimentID;             //!< e_RDCSourceExperimentID;
      const StarTag e_RDCAuthAsymID1;                    //!< e_RDCAuthAsymID1;
      const StarTag e_RDCAuthSeqID1;                     //!< e_RDCAuthSeqID1;
      const StarTag e_RDCAuthCompID1;                    //!< e_RDCAuthCompID1;
      const StarTag e_RDCAuthAtomID1;                    //!< e_RDCAuthAtomID1;
      const StarTag e_RDCAuthAsymID2;                    //!< e_RDCAuthAsymID2;
      const StarTag e_RDCAuthSeqID2;                     //!< e_RDCAuthSeqID2;
      const StarTag e_RDCAuthCompID2;                    //!< e_RDCAuthCompID2;
      const StarTag e_RDCAuthAtomID2;                    //!< e_RDCAuthAtomID2;
      const StarTag e_RDCEntryID;                        //!< e_RDCEntryID;
      const StarTag e_RDCRDCConstraintListID;            //!< e_RDCRDCConstraintListID;

      //! RDC list tags
      const StarTag e_RDCListSfCategory;                 //!< e_RDCListSfCategory;
      const StarTag e_RDCListEntryID;                    //!< e_RDCListEntryID;
      const StarTag e_RDCListID;                         //!< e_RDCListID;
      const StarTag e_RDCListConstraintFileID;           //!< e_RDCListConstraintFileID;
      const StarTag e_RDCListBlockID;                    //!< e_RDCListBlockID;
      const StarTag e_RDCListDetais;                     //!< e_RDCListDetais;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StarTags();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief add a new tag
      //! @param DESCRIPTOR name of tag
      //! @param TAG_DATA data for the tag
      //! @return the created enum for that tag
      StarTag AddStarTag( const std::string &DESCRIPTOR, const StarTagData &TAG_DATA);

      //! @brief get a tag from a string
      //! @param CATEGORY this tag belongs to
      //! @param DESCRIPTION string description for the enum
      //! @return a tag from a string
      static StarTag GetTagFromString( const StarTagCategory &CATEGORY, const std::string &DESCRIPTION);

      //! @brief read Star file from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      static storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > ReadStarFile( std::istream &ISTREAM);

      //! @brief gets a chain id from a string since the star format allows authors to identify chains however they want
      //! @param CURRENT_ID current string to be converted
      //! @param PREVIOUS_IDS previously seen strings mapped to chain ids
      //! @return chain id
      static char GetChainID( const std::string &CURRENT_ID, storage::Map< std::string, char> &PREVIOUS_IDS);

    }; // class StarTags

    //! @brief get access to all star tags
    BCL_API const StarTags &GetStarTags();

  } // namespace nmr

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< nmr::StarTagData, nmr::StarTags>;

  } // namespace util
} // namespace bcl

#endif // BCL_NMR_STAR_TAGS_H_
