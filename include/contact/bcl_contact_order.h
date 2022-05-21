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

#ifndef BCL_CONTACT_ORDER_H_
#define BCL_CONTACT_ORDER_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "score/bcl_score_protein_model.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Order
    //! @brief This is a Function derived class for calculating the contact order of a ProteinModel
    //! @details Contact order represents the average sequence separation between physically contacting residue pairs
    //! in a given amino acid sequence. The higher value indicates the existence of long distance contacts while a
    //! a small value indicates mostly local short distance contacts.
    //!
    //! @see @link example_contact_order.cpp @endlink
    //! @author woetzen, karakam
    //! @date 12.07.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Order :
      public score::ProteinModel
    {
    public:

    ///////////
    // types //
    ///////////

      //! @brief type of normalization that can be applied to contact order
      enum NormalizationType
      {
        e_Absolute,
        e_RelativeAAsUsed,
        e_RelativeSequenceLength,
        e_RelativeSqrAAsUsed,
        e_RelativeSqrSequenceLength,
        s_NumberNormalizationType
      };

      //! @brief NormalizationType as string
      //! @param NORMALIZATION_TYPE the NormalizationType
      //! @return the string for the NormalizationType
      static const std::string &GetNormalizationTypeDescriptor( const NormalizationType &NORMALIZATION_TYPE);

      //! @brief NormalizationTypeEnum is used for I/O of NormalizationType
      typedef
        util::WrapperEnum
        <
          NormalizationType,
          &GetNormalizationTypeDescriptor,
          s_NumberNormalizationType
        > NormalizationTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! type of normalization
      NormalizationTypeEnum m_Normalization;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! bool: true to cache neighbor lists. This speeds up calculations when contact recovery and other scores that
      //! depend on the neighbor lists are also being computed
      bool m_CacheNeighborLists;

      //! generator to be used
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > m_NeighborGenerator;

      static math::FunctionCached< assemble::ProteinModel, double> s_AbsoluteInstanceCached;

    public:

    //////////
    // data //
    //////////

      //! sequence distance of two aas to be considered for contactorder
      static const size_t s_ContactOrderSequencesDistanceCutoff = 4;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Order();

      //! @brief constructor from normalization type, scheme string and whether to use cache generating neighbor lists
      //! @brief NORMALIZE type of normalization
      //! @brief SCHEME scheme to report
      //! @brief CACHE whether to use cache generating neighbor lists
      Order( const NormalizationType NORMALIZE, const std::string &SCHEME, const bool CACHE = false);

      //! @brief copy constructor
      Order *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief access to the normalization type used
      //! @return the normalization enum used
      const NormalizationTypeEnum &GetNormalization() const
      {
        return m_Normalization;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief contact order for aa sequence
      //! @param AA_SEQUENCE sequence of interest
      //! @return contact order of sequence
      double ContactOrder( const biol::AASequence &AA_SEQEUNCE) const;

      //! @brief contact order for the given chain
      //! @param CHAIN chain of interest
      //! @return contact order for given chain
      double ContactOrder( const assemble::Chain &CHAIN) const;

      //! @brief get score type
      //! @return score type
      score::ProteinModel::Type GetType() const
      {
        return score::ProteinModel::e_Structure;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the contact order for a given protein model
      //! @param PROTEIN_MODEL protein model
      //! @return contact order for the given protein model
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    private:

    ////////////
    // helper //
    ////////////

      //! @brief normalize the contact order
      //! @param CONTACT_ORDER the unnormalized contact order
      //! @param NR_AAS_USED nr of amino acids used
      //! @param SEQUENCE_LENGTH the length of the amino acid sequence
      //! @return normalized contact order, depending of m_Normalization
      double Normalize( const double CONTACT_ORDER, const size_t NR_AAS_USED, const size_t SEQUENCE_LENGTH) const;

      //! @brief calculate contact order form NeighborListContainer
      //! @param NEIGHBOR_LIST_CONTAINER container of neighbors
      //! @return the contact order
      static double ContactOrder( const assemble::AANeighborListContainer &NEIGHBOR_LIST_CONTAINER);

    }; // class Order

  } // namespace contact
} // namespace bcl

#endif //BCL_CONTACT_ORDER_H_
