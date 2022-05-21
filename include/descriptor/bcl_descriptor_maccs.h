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

#ifndef BCL_DESCRIPTOR_MACCS_H_
#define BCL_DESCRIPTOR_MACCS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "graph/bcl_graph_const_graph.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MACCS
    //! @brief MACCS key fingerprint, based off PubChem Maccs keys; 938 bits long
    //!
    //! @see @link example_descriptor_maccs.cpp @endlink
    //! @see @link ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt @endlink
    //! @author perrye, mendenjl, geanesar
    //! @date Jun 25, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MACCS :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {
    public:

      //! Number of elements and smart queries considered by MACCS keys
      enum
      {
        s_NumberElements = 92,
        s_NumberSMARTSQueries = 554
      };

    private:

      static const size_t      s_ElementArray[ s_NumberElements];
      static const std::string s_AroSMARTSArray[ s_NumberSMARTSQueries];

      storage::Vector< size_t>                    m_ElementHistogram; //!< Histogram of element type
      storage::Vector< size_t>                    m_RingSizeHistogram; //!< Histogram of ring sizes

      //! Histogram of ring sizes of types pure carbon, carbon + nitrogen, and heteroatoms
      storage::Vector< storage::Vector< size_t> > m_RingTypeHistogram;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MACCS();

      //! @brief virtual copy constructor
      MACCS *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 938;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class MACCS

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MACCS_H_
