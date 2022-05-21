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

#ifndef BCL_DESCRIPTOR_AA_SS_TM_PREDICTION_H_
#define BCL_DESCRIPTOR_AA_SS_TM_PREDICTION_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASSTMPrediction
    //! @brief generates state secondary structure prediction for a given amino acid
    //! @details This Interface derived class allows generation of descriptions from 9 state secondary
    //! structure/transmembrane prediction of a given amino acid. It concatenates these descriptions for each of the
    //! given SSMethods and returns them.
    //!
    //! @see @link example_descriptor_aa_ss_tm_prediction.cpp @endlink
    //! @author karakam, mendenjl
    //! @date Apr 24, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASSTMPrediction :
      public BaseElement< biol::AABase, float>
    {

    private:

    //////////
    // data //
    //////////

      //! SSMethod that will be used
      sspred::Method m_Method;

      //! if true, consider membrane states (gives 9 states instead of 3)
      bool m_ConsiderMembraneStates;

      //! if true, consider secondary structure states
      bool m_ConsiderSecondaryStructure;

      //! min sse sizes for consideration
      size_t m_MinHelixSize;      //!< Minimum number of residues to consider an SSE a helix, if smaller, return coil
      size_t m_MinStrandSize;     //!< Minimum number of residues to consider an SSE a strand; if smaller, return coil
      size_t m_MinMembraneSpan;   //!< Minimum number of residues in the membrane; if smaller, return solution

      //! @brief alias, which is dependent on all other member variables
      std::string m_Alias;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from an alias and what states to consider
      //! @param ALIAS alias to use over the command line for this class
      //! @param CONSIDER_MEMBRANE whether to consider the membrane states
      //! @param CONSIDER_SS whether to consider secondary structure states
      AASSTMPrediction
      (
        const std::string &ALIAS,
        const bool CONSIDER_MEMBRANE,
        const bool CONSIDER_SS
      );

      //! @brief Clone function
      //! @return pointer to new AASSTMPrediction
      AASSTMPrediction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief set the method used for calculation
      //! @param METHOD the new method to use
      void SetMethod( const sspred::Method &METHOD);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AASSTMPrediction

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_SS_TM_PREDICTION_H_
