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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_MOMENT_OF_INERTIA_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_MOMENT_OF_INERTIA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_aa_exposure_interface.h"
#include "biol/bcl_biol_aa_type_data.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelMomentOfInertia
    //! @brief calculate the moments of inertia and the moment of inertia tensor for a protein model
    //!
    //! @see @link example_assemble_protein_model_moment_of_inertia.cpp @endlink
    //! @author woetzen
    //! @date Apr 16, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelMomentOfInertia :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! aatype property used as weight
      biol::AATypeData::PropertyTypeEnum m_PropertyType;

      //! exposure weight function
      util::ShPtr< AAExposureInterface> m_ExposureWeight;

      //! consider only negative transfer free energies
      bool m_ConsiderOnlyNegativeEnergies;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param PROPERTY the property to use as weight
      //! @param EXPOSURE_WEIGHT_FUNCTION multiply each property with aa exposure
      //! @param CONSIDER_ONLY_NEGATIVE_PROPERTY consider only negative energies in calculations
      ProteinModelMomentOfInertia
      (
        const biol::AATypeData::PropertyType &PROPERTY,
        const util::ShPtr< AAExposureInterface> &EXPOSURE_WEIGHT_FUNCTION,
        const bool CONSIDER_ONLY_NEGATIVE_PROPERTY
      );

      //! @brief Clone function
      //! @return pointer to new ProteinModelMomentOfInertia
      ProteinModelMomentOfInertia *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief return a 4*n matrix with 3 coordinates and exposure * weight in each row (4 columns) and number amino acid rows
      //! @param PROTEIN_MODEL protein model
      //! @return 4 * n matrix
      linal::Matrix< double> ProteinModelAAExposureToCoordinateWeightMatrix( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief return a 4*n matrix with 3 coordinates and weight in each row (4 columns) and number amino acid rows
      //! @param PROTEIN_MODEL protein model
      //! @return 4 * n matrix
      linal::Matrix< double> ProteinModelToCoordinateWeightMatrix( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief calculate transformation that translates into the center of weights and that sorts principal axes of inertia according to principal moments of inertia x - smallest, z - largest
      //! @param PROTEIN_MODEL protein model
      //! @return transformation matrix and principal moments of inertias
      storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
      TransformationAndMoments( const ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ProteinModelMomentOfInertia

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_MODEL_MOMENT_OF_INERTIA_H_
