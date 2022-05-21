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

#ifndef BCL_DESCRIPTOR_AA_INFO_FROM_SYM_MATRIX_FILE_H_
#define BCL_DESCRIPTOR_AA_INFO_FROM_SYM_MATRIX_FILE_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_pair.h"
#include "bcl_descriptor_type.h"
#include "storage/bcl_storage_symmetric_matrix.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAInfoFromSymMatrixFile
    //! @brief AAInfoFromSymMatrixFiles the output of descriptors into one vector
    //!
    //! @see @link example_descriptor_aa_info_from_sym_matrix_file.cpp @endlink
    //! @author teixeipl, mendenjl
    //! @date Feb 01, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAInfoFromSymMatrixFile :
      public BasePair< biol::AABase, float>
    {
    private:

    //////////
    // data //
    //////////

      //! Symmetric matrix storing data from file
      storage::SymmetricMatrix< double> m_SymmetricMatrix;

      //! Suffix for symmetric matrix BCL file that contains data
      std::string m_FileSuffix;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      AAInfoFromSymMatrixFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
      //! @param STORAGE storage for the descriptor
      //! @return true, if the calculation was successful
      virtual void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT_A,
        const iterate::Generic< const biol::AABase> &ELEMENT_B,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    }; // class AAInfoFromSymMatrixFile

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_INFO_FROM_SYM_MATRIX_FILE_H_
