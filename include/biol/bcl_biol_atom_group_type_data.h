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

#ifndef BCL_BIOL_ATOM_GROUP_TYPE_DATA_H_
#define BCL_BIOL_ATOM_GROUP_TYPE_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "restraint/bcl_restraint_sas_data_parameters.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomGroupTypeData
    //! @brief This is a low level helper class to store biological atom group properties
    //! @details This class acts as the storage class for the enumerator AtomGroups
    //!
    //! @see @link example_biol_atom_group_type_data.cpp @endlink
    //! @author putnamdk
    //! @date 11/21/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomGroupTypeData :
      public math::FunctionInterfaceSerializable< restraint::SasDataParameters, double>

      {
      friend class AtomGroupTypes;

    private:

    //////////
    // data //
    //////////

      //! effective form factor for in vacuo
      util::Implementation< math::FunctionInterfaceSerializable< restraint::SasDataParameters, double> > m_Vacuo;

      //! effective form factor for water
      util::Implementation< math::FunctionInterfaceSerializable< restraint::SasDataParameters, double> > m_Water;

      std::string m_BaseAtomType;
      double m_DisplacedSolventVolume;
      double m_Radius;
      size_t m_HydrogenCount;
      double m_H2OScatteringLength;
      double m_D2OScatteringLength;

      //! effective mass of the atom group
      double m_Mass;

      //! effective surface area
      double m_SurfaceArea;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! solvent density
      static const double s_SolventDensity;

      //////////////////////////////////
      // construction and destruction //group
      //////////////////////////////////

      //! @brief construct undefined atom group type
      AtomGroupTypeData();

      //! @brief construct atom type from given information
      //! @param BASE_ATOM_TYPE base atom type
      //! @param DISPLACED_SOLVENT_VOLUME Volume of solvent displaced by atom group
      //! @param RADIUS Radius of theoretical sphere of displaced solvent
      //! @param HYDROGEN_COUNT Number of Hydrogens bound to Heavy Atom
      AtomGroupTypeData
      (
        const AtomType &BASE_ATOM_TYPE,
        const double &DISPLACED_SOLVENT_VOLUME,
        const double &RADIUS,
        const size_t &HYDROGEN_COUNT,
        const double &H2O_SCATTERING_LENGTH,
        const double &D20_SCATTERING_LENGTH
      );

      //! @brief virtual copy constructor
      AtomGroupTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator
      //! @param FORM_FACTOR_DATA Data to Calculate the form factors
      //! @return form factor
      double operator()( const restraint::SasDataParameters &FORM_FACTOR_DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class AtomGroupTypeData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ATOM_GROUP_TYPE_DATA_H_

