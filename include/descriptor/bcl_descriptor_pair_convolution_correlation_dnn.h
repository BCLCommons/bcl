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

#ifndef BCL_DESCRIPTOR_PAIR_CONVOLUTION_CORRELATION_DNN_H_
#define BCL_DESCRIPTOR_PAIR_CONVOLUTION_CORRELATION_DNN_H_

// include the namespace header
#include "bcl_descriptor.h"
#include "bcl_descriptor_combine.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PairConvolutionCorrelationDNN
    //! @brief Estimates binding affinity (pKd units) in a pose-dependent manner with a deep neural network
    //!
    //! @see @link example_descriptor_pair_convolution_correlation_dn_n.cpp @endlink
    //! @author brownbp1
    //! @date May 19, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PairConvolutionCorrelationDNN :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    public:

    //////////
    // data //
    //////////

      // enables weighting of pKd by applicability domain score
      bool m_ApplicabilityWeight;

      // molecule of interest for binding affinity calculation
      chemistry::FragmentComplete m_Molecule;

      // SDF MDL property name for pocket of interest
      std::string m_MDLPocketName;

      // filename for protein binding pocket
      std::string m_ProteinPocketFilename;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PairConvolutionCorrelationDNN();

      //! @brief specify applicability domain model constructor
      explicit PairConvolutionCorrelationDNN
      (
        const bool &APPLICABILITY_WEIGHT = false
      );

      //! @brief detailed constructor
      explicit PairConvolutionCorrelationDNN
      (
        const std::string &MDL_PROPERTY = std::string(),
        const std::string &POCKET_FILENAME = std::string(),
        const bool &APPLICABILITY_WEIGHT = false
      );

      //! @brief Clone function
      //! @return pointer to new PairConvolutionCorrelationDNN
      PairConvolutionCorrelationDNN *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class PairConvolutionCorrelationDNN

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_PAIR_CONVOLUTION_CORRELATION_DNN_H_
