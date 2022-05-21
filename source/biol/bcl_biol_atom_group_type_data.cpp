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
#include "biol/bcl_biol_atom_group_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_sum_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomGroupTypeData::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomGroupTypeData())
    );

    //! solvent density
    const double AtomGroupTypeData::s_SolventDensity( 0.334);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined atom type
    AtomGroupTypeData::AtomGroupTypeData() :
      m_Vacuo(),
      m_Mass( util::GetUndefined< double>()),
      m_SurfaceArea( util::GetUndefined< double>())
    {
    }

    //! @brief construct atom type from given information
    //! @param BASE_ATOM_TYPE base atom type
    //! @param DISPLACED_SOLVENT_VOLUME Volume of solvent displaced by atom group
    //! @param RADIUS Radius of theoretical sphere of displaced solvent
    //! @param HYDROGEN_COUNT Number of Hydrogens bound to Heavy Atom
    AtomGroupTypeData::AtomGroupTypeData
    (
      const AtomType &BASE_ATOM_TYPE,
      const double &DISPLACED_SOLVENT_VOLUME,
      const double &RADIUS,
      const size_t &HYDROGEN_COUNT,
      const double &H2O_SCATTERING_LENGTH,
      const double &D2O_SCATTERING_LENGTH
    ) :
      m_Vacuo(),
      m_Water(),
      m_BaseAtomType( BASE_ATOM_TYPE.GetName()),
      m_DisplacedSolventVolume( DISPLACED_SOLVENT_VOLUME),
      m_Radius( RADIUS),
      m_HydrogenCount( HYDROGEN_COUNT),
      m_H2OScatteringLength( H2O_SCATTERING_LENGTH),
      m_D2OScatteringLength( D2O_SCATTERING_LENGTH),
      m_Mass( DISPLACED_SOLVENT_VOLUME * s_SolventDensity),
      m_SurfaceArea( ( -math::Pow( DISPLACED_SOLVENT_VOLUME, 2.0 / 3.0)) / ( 4.0 * math::g_Pi))
    {

      // initialize structure factors for water
      util::ShPtr< math::SumFunction< restraint::SasDataParameters, double> > water_factors
      (
        new math::SumFunction< restraint::SasDataParameters, double>()
      );

      // get the crommer mann constants for hydrogen
      static const math::SumFunction< restraint::SasDataParameters, double>
      h_form_factor( GetAtomTypes().H->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // get the crommer mann constants for oxygen
      static const math::SumFunction< restraint::SasDataParameters, double>
      o_form_factor( GetAtomTypes().O->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // add hydrogen structure factors to water_factors
      *water_factors += double( 2.0) * h_form_factor;

      // add oxygen structure factors to water_factors
      *water_factors += double( 1.0) * o_form_factor;
      m_Water = *water_factors;

      // initialize structure factors
      math::SumFunction< restraint::SasDataParameters, double> structure_factors;

      // initialize with element structure factor
      structure_factors += BASE_ATOM_TYPE->GetElementType()->GetStructureFactor();

      // add hydrogen structure factors (use h_form_factor from above)
      structure_factors += double( HYDROGEN_COUNT) * h_form_factor;

      // set vacuo form factor
      m_Vacuo = structure_factors;
    }

    //! @brief virtual copy constructor
    AtomGroupTypeData *AtomGroupTypeData::Clone() const
    {
      return new AtomGroupTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomGroupTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator
    //! @param FORM_FACTOR_DATA Data to Calculate the form factors
    //! @return form factor
    double AtomGroupTypeData::operator()( const restraint::SasDataParameters &FORM_FACTOR_DATA) const
    {

      const double &Qvalue( FORM_FACTOR_DATA.GetQValue());
      const double &Sasa( FORM_FACTOR_DATA.GetSasaValue() / 100);
      const double &ExcludedVolumeParameter( FORM_FACTOR_DATA.GetExcludedVolume());
      const double &HydrationShellParameter( FORM_FACTOR_DATA.GetHydrationShell());
      const double &DeteriumPercentage( FORM_FACTOR_DATA.GetDeuteriumExchangeRate());

      double FormFactor( util::GetUndefinedDouble());

      if( FORM_FACTOR_DATA.GetSansImplementation())
      {
        double water_contribution( -0.00562);
        double deuterium_contribution( 0.0697 * FORM_FACTOR_DATA.GetDeuteriumExchangeRate());
        double solvent_scattering_length_density( water_contribution + deuterium_contribution);

        double vacuumScatteringLength( 0.0);

        // Hydrodens bound to Carbon will not exchange
        if( m_BaseAtomType == "C")
        {
          vacuumScatteringLength = m_H2OScatteringLength;
        }

        // Perform the deuterium exchange for Nitrogen, Oxygen, and Sulfer
        if( m_BaseAtomType == "N" || m_BaseAtomType == "O" || m_BaseAtomType == "S")
        {
          vacuumScatteringLength = ( m_H2OScatteringLength * ( 1 - DeteriumPercentage)) + ( m_D2OScatteringLength * DeteriumPercentage);
        }

        if( m_HydrogenCount == 0)
        {
          vacuumScatteringLength = m_H2OScatteringLength;
        }

        // Perform Calculation of water factor based on Deterium Exchange rate

        double water( -0.168);
        double deuterium( 1.915);

        double hydration_scattering_length( water * ( 1 - DeteriumPercentage) + ( deuterium * DeteriumPercentage));

        FormFactor =
          vacuumScatteringLength -
          ExcludedVolumeParameter *solvent_scattering_length_density * m_DisplacedSolventVolume * std::exp( math::Sqr( Qvalue) * m_SurfaceArea) +
          HydrationShellParameter *Sasa *hydration_scattering_length;
      }
      else
      {
        FormFactor =
          m_Vacuo->operator ()( FORM_FACTOR_DATA) -
          ExcludedVolumeParameter *m_Mass * std::exp( math::Sqr( Qvalue) * m_SurfaceArea) +
          HydrationShellParameter *Sasa *m_Water->operator()( FORM_FACTOR_DATA);
      }

      return FormFactor;
    }
    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomGroupTypeData::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Stores biological and group properties of atom groups.");
      serializer.AddInitializer
      (
        "base atom type",
        "",
        io::Serialization::GetAgent( &m_BaseAtomType)
      );
      serializer.AddInitializer
      (
        "displaced solvent volume",
        "",
        io::Serialization::GetAgent( &m_DisplacedSolventVolume)
      );
      serializer.AddInitializer
      (
        "radius",
        "",
        io::Serialization::GetAgent( &m_Radius)
      );
      serializer.AddInitializer
      (
        "hydrogen count",
        "",
        io::Serialization::GetAgent( &m_HydrogenCount)
      );
      serializer.AddInitializer
      (
        "h2o scattering length",
        "",
        io::Serialization::GetAgent( &m_H2OScatteringLength)
      );
      serializer.AddInitializer
      (
        "s2o scattering length",
        "",
        io::Serialization::GetAgent( &m_D2OScatteringLength)
      );
      serializer.AddDataMember
      (
        "mass",
        io::Serialization::GetAgent( &m_Mass)
      );
      serializer.AddDataMember
      (
        "surface area",
        io::Serialization::GetAgent( &m_SurfaceArea)
      );
      serializer.AddDataMember
      (
        "form factor vacuo",
        io::Serialization::GetAgent( &m_Vacuo)
      );
      serializer.AddDataMember
      (
        "form factor water",
        io::Serialization::GetAgent( &m_Water)
      );

      return serializer;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomGroupTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Vacuo, ISTREAM);
      io::Serialize::Read( m_Mass, ISTREAM);
      io::Serialize::Read( m_SurfaceArea, ISTREAM);
      io::Serialize::Read( m_Water, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AtomGroupTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Vacuo, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Mass, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SurfaceArea, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Water, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
