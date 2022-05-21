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
#include "descriptor/bcl_descriptor_aa_info_from_sym_matrix_file.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAInfoFromSymMatrixFile::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAInfoFromSymMatrixFile()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAInfoFromSymMatrixFile *AAInfoFromSymMatrixFile::Clone() const
    {
      return new AAInfoFromSymMatrixFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAInfoFromSymMatrixFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAInfoFromSymMatrixFile::GetAlias() const
    {
      static const std::string s_name( "AAInfoFromSymMatrixFile");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAInfoFromSymMatrixFile::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
    //! @param STORAGE storage for the descriptor
    void AAInfoFromSymMatrixFile::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT_A,
      const iterate::Generic< const biol::AABase> &ELEMENT_B,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( size_t( std::max( ELEMENT_A->GetSeqID(), ELEMENT_B->GetSeqID())) > m_SymmetricMatrix.GetSize())
      {
        STORAGE = util::GetUndefined< float>();
      }
      else
      {
        STORAGE( 0) = m_SymmetricMatrix( ELEMENT_A->GetSeqID() - 1, ELEMENT_B->GetSeqID() - 1);
      }
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AAInfoFromSymMatrixFile::SetObjectHook()
    {
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // Get the filename
      util::ShPtr< util::Wrapper< std::string> > sp_filename_wrapper
      (
        sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
      );

      std::string pdb_filename( sp_filename_wrapper->GetData());

      // Remove the last extension
      std::string basename( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( pdb_filename)));

      // Add my extension
      std::string sym_matrix_filename( basename + m_FileSuffix);
      // Read in the file's data
      io::IFStream input;
      io::File::MustOpenIFStream( input, sym_matrix_filename);
      // Read and close file stream
      io::Serialize::Read( m_SymmetricMatrix, input);
      io::File::CloseClearFStream( input);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAInfoFromSymMatrixFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "data for amino acid pairs based on their sequence index and the given symmetric matrix file"
      );
      parameters.AddInitializer
      (
        "suffix", "suffix of the bcl symmetric matrix file to be read in",
        io::Serialization::GetAgent( &m_FileSuffix)
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
