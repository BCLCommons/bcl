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

#ifndef BCL_NMR_SPECTRUM_H_
#define BCL_NMR_SPECTRUM_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_nmr_signal.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Spectrum
    //! @brief This class is used to store spectra (ensembles of signals). The signals consist of a
    //! shared pointer vector of signals1D with the same dimension as the spectrum. Several
    //! are combined to a spectrum. A molecule can have multiple spectra.
    //!
    //! @see @link example_nmr_spectrum.cpp @endlink
    //! @author mueller, butkiem1
    //! @date 04/11/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Spectrum :
      public util::ObjectInterface
    {

    public:

      //! @brief enumerator for all NMR spectra types
      enum SpecType
      {
        e_NoSpecType,
        e_1H,
        e_11B,
        e_13C,
        e_15N,
        e_31P,
        e_COSY,
        e_DEPT,
        e_HETCOR,
        e_NOESY,
        s_NumberSpecTypes
      };

      //! @brief SpecType as string
      //! @param SPEC_TYPE the SpecType
      //! @return the string for the SpecType
      static const std::string &GetSpecTypeDescriptor( const SpecType &SPEC_TYPE);

      //! @brief SpecTypeEnum is used for I/O of SpecType
      typedef util::WrapperEnum< Spectrum::SpecType, &Spectrum::GetSpecTypeDescriptor, Spectrum::s_NumberSpecTypes> SpecTypeEnum;

    private:

    //////////
    // data //
    //////////

      SpecTypeEnum               m_Type;             //!< type of spectrum
      util::ShPtrVector< Signal> m_Signals;          //!< list of signals
      double                     m_Temperature;      //!< temperature of spectrum
      std::string                m_Solvent;          //!< solvent of spectrum
      double                     m_FieldStrength;    //!< field strength of spectrum in MHz
      std::string                m_AssignmentMethod; //!< method used for spectra assignment, e.g. "Calculated", "Manual", "Automatic"

    public:

      static const std::string s_FieldStrengthPropertyDescriptor;
      static const std::string s_TemperaturePropertyDescriptor;
      static const std::string s_SolventPropertyDescriptor;
      static const std::string s_AssignmentMethodDescriptor;
      static const std::string s_UnreportedDescriptor;
      static const std::string s_SpectrumIdentifier;
      static const char        s_SpectraDescriptorSeparator = ':';
      static const char        s_SignalSeparator = '|';
      static const char        s_PropertySeparator = ' ';

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief standard constructor
      Spectrum();

      //! @brief construct spectrum from its properties
      Spectrum( const SpecTypeEnum &TYPE, const util::ShPtrVector< Signal> &SIGNALS);

      //! @brief virtual copy constructor
      Spectrum *Clone() const
      {
        return new Spectrum( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief  get type of spectrum
      //! @return type of spectrum
      const SpecTypeEnum &GetSpecType() const
      {
        return m_Type;
      }

      //! @brief set type of spectrum to given TYPE
      void SetType( const SpecTypeEnum &TYPE)
      {
        m_Type = TYPE;
      }

      //! @brief get Signals
      //! @return Signals
      const util::ShPtrVector< Signal> &GetSignals() const
      {
        return m_Signals;
      }

      //! @brief Set signals to given SIGNALS
      void SetSignals( const util::ShPtrVector< Signal> &SIGNALS)
      {
        m_Signals = SIGNALS;
      }

      //! @brief get temperature
      //! @return temperature
      const double &GetTemperature() const
      {
        return m_Temperature;
      }

      //! @brief set temperature to given TEMPERATURE
      //! @param TEMPERATURE temperature
      void SetTemperature( const double &TEMPERATURE)
      {
        m_Temperature = TEMPERATURE;
      }

      //! @brief get field strength
      //! @return field strength
      const double &GetFieldStrength() const
      {
        return m_FieldStrength;
      }

      //! @brief set field strength to given FIELDSTRENGHT
      //! @param FIELDSTRENGTH field strength
      void SetFieldStrength( const double &FIELDSTRENGTH)
      {
        m_FieldStrength = FIELDSTRENGTH;
      }

      //! @brief get solvent
      //! @return solvent as a string
      const std::string &GetSolvent() const
      {
        return m_Solvent;
      }

      //! @brief set assignment method to given
      //! @param METHOD assignment method
      void SetAssignmentMethod( const std::string &METHOD)
      {
        m_AssignmentMethod = METHOD;
      }

      //! @brief get assignment method
      //! @return assignment method as a string
      const std::string &GetAssignmentMethod() const
      {
        return m_AssignmentMethod;
      }

      //! @brief set solvent to given solvent
      //! @param SOLVENT solvent
      void SetSolvent( const std::string &SOLVENT)
      {
        m_Solvent = SOLVENT;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief method for adding a spectrum to a small molecule
      //! @param MOLECULE the small molecule the spectrum is added to
      void AddSpectrum( chemistry::ConformationInterface &MOLECULE) const;

      //! @brief get all signals involving given atom
      //! @param ATOM atom of interest
      //! @return util::SiPtrVector< const Signal> all signals involving atom of interest
      util::SiPtrVector< const Signal> DetermineInvolvedSignals( const chemistry::AtomConformationalInterface &ATOM) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read Spectrum from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write Spectrum to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief convert property line to property map
      //! @param CONFORMATION the conformation that should have the given property
      //! @param PROPERTY_LINE string of the property line "0:Unreported 1:CDCl3" create a map of size_t:string
      //! @return map that maps id to property
      static storage::Map< size_t, std::string> GetPropertyMap
      (
        const chemistry::ConformationInterface &CONFORMATION,
        const std::string &PROPERTY_LINE
      );

      //! @brief convert property map to property line "0:Unreported 1:CDCl3"
      //! @param PROPERTY_MAP that maps size_t to string
      //! @return string of the form "0:Unreported 1:CDCl3"
      static std::string MapToLine( const storage::Map< size_t, std::string> &PROPERTY_MAP);

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief generate the spectra
      //! @param MOLECULE the small molecule associated with the spectra
      //! @return a Map of spectra, where the key is the spectra number
      static storage::Map< size_t, util::ShPtr< Spectrum> > GenerateSpectra
      (
        const chemistry::ConformationInterface &MOLECULE
      );

      //! @brief read signals from and MDL Spectrum line
      //! @param MOLECULE the small molecule the spectrum is for
      //! @param SPECTRUM_LINE the spectrum line with all the signals
      //! @return a ShPtrVector of Signals
      static util::ShPtrVector< Signal> ReadSignals
      (
        const chemistry::ConformationInterface &MOLECULE,
        const std::string &SPECTRUM_LINE
      );

    }; // class Spectrum

  } // namespace nmr
} // namespace bcl

#endif //BCL_NMR_SPECTRUM_H_

