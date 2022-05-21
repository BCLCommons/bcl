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

#ifndef BCL_DESCRIPTOR_MOLECULE_ENTROPY_QHA_H_
#define BCL_DESCRIPTOR_MOLECULE_ENTROPY_QHA_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeEntropyQHA
    //! @brief Quasi-harmonic approximation of the conformational entropy of a small molecule
    //! @details Relates two conformational ensembles (i.e. a "local" ensemble and a "global" ensemble) through
    //!          principal component analysis (PCA) of the coordinate space. Unless explicitly provided, the ensembles
    //!          are obtained by sampling conformational space with chemistry::SampleConformations. Local sampling
    //!          of a small molecule 3D conformer occurs through uniform and Gaussian dihedral angle perturbations.
    //!          Local sampling occurs almost entirely within the starting dihedral bin for each considered rotamer
    //!          (2 standard deviations, 95%). Global sampling utilizes a full conformational search of the molecule.
    //!          Descriptor output indices correspond to the global_s, local_s, global_s - local_s, and -ln(local_s/global_s).
    //!          The first 4 indices are actual entropy estimates, while the latter 4 indices are the PCA eigenvalue sums.
    //!
    //!          The QHA equation we referenced is eq. 13 in "On the calculation of entropy from covariance matrices
    //!          of the atomic fluctuations", Andriciaei and Karplus, 2001, The Journal of Chemical Physics.
    //!          https://doi.org/10.1063/1.1401821
    //!
    //! @see @link example_descriptor_molecule_entropy_qha.cpp @endlink
    //! @author brownbp1
    //! @date Aug 22, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeEntropyQHA :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {
    public:

    private:
      // Sample molecule conformers
      chemistry::SampleConformations m_SampleConfs;

      // Local conformational ensemble filename
      std::string m_LocalEnsembleFilename;

      // Global conformational ensemble filename
      std::string m_GlobalEnsembleFilename;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeEntropyQHA();

      //! @brief construct with conformation sampler options
      explicit MoleculeEntropyQHA
      (
        const chemistry::SampleConformations &SAMPLER
      );

      //! @brief copy constructor
      //! @return a pointer to a copy of this class
      MoleculeEntropyQHA *Clone() const;

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
        return size_t( 8);
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

    }; // class MoleculeEntropyQHA

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_ENTROPY_QHA_H_
