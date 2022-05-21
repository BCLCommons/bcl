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
#include "nmr/bcl_nmr_spectrum.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief SpecType as string
    //! @param SPEC_TYPE the SpecType
    //! @return the string for the SpecType
    const std::string &Spectrum::GetSpecTypeDescriptor( const SpecType &SPEC_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "NO_SPECTYPE",
        "1H",
        "11B",
        "13C",
        "15N",
        "31P",
        "COSY",
        "DEPT",
        "HETCOR",
        "NOESY",
        GetStaticClassName< SpecType>()
      };

      return s_descriptors[ SPEC_TYPE];
    }

  //////////
  // data //
  //////////

    const std::string Spectrum::s_FieldStrengthPropertyDescriptor( "Field Strength [MHz]");
    const std::string Spectrum::s_TemperaturePropertyDescriptor( "Temperature [K]");
    const std::string Spectrum::s_SolventPropertyDescriptor( "Solvent");
    const std::string Spectrum::s_AssignmentMethodDescriptor( "Assignment Method");
    const std::string Spectrum::s_UnreportedDescriptor( "Unreported");
    const std::string Spectrum::s_SpectrumIdentifier( "Spectrum");

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Spectrum::s_Instance
    (
      GetObjectInstances().AddInstance( new Spectrum())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    Spectrum::Spectrum()
    {
    }

    //! @brief construct spectrum from its properties
    Spectrum::Spectrum( const SpecTypeEnum &TYPE, const util::ShPtrVector< Signal> &SIGNALS) :
      m_Type( TYPE),
      m_Signals( SIGNALS),
      m_Temperature( 298),
      m_Solvent(),
      m_FieldStrength( util::GetUndefined< double>())
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief method for adding a spectrum to a small molecule
    //! @param MOLECULE the small molecule the spectrum is added to
    void Spectrum::AddSpectrum( chemistry::ConformationInterface &MOLECULE) const
    {
      std::string value;

      // for every signal
      for
      (
        util::ShPtrVector< Signal>::const_iterator
          itr_signal( m_Signals.Begin()), itr_signal_end( m_Signals.End());
        itr_signal != itr_signal_end;
        ++itr_signal
      )
      {
        const Signal1D &current_signal( *( *itr_signal)->GetSignals1D().FirstElement());

        // find the position of the atom that belongs to the given reference in signal1D object
        const size_t atom_index( MOLECULE.GetAtomIndex( *current_signal.GetAtomInvolvedInSignal()));

        const storage::Pair< Signal1D::PatternTypeEnum, double> &current_pattern( current_signal.GetPattern()( 0));

        // write spectrum    structure :  shift;couplingconstPattern;AtomIndex|
        value += util::Format().FFP( 1)( current_signal.GetChemicalShift()); // shift
        value += Signal1D::s_InformationSeparator;                   // separator
        value += util::Format().FFP( 1)( current_pattern.Second());  // coupling constant
        value += Signal1D::GetPatternTypeSingleLetterCode( current_pattern.First()); // coupling pattern
        value += Signal1D::s_InformationSeparator;                   // separator
        value += util::Format()( atom_index);                        // index of atom
        value += s_SignalSeparator;
      }

      // construct spectrum name from type and count
      const std::string spectrum_name_base( s_SpectrumIdentifier + " " + m_Type.GetString() + " ");

      // increase count until the spectrum does not exist as a misc property
      const storage::Map< size_t, util::ShPtr< Spectrum> > spectra( GenerateSpectra( MOLECULE));
      const size_t count( spectra.IsEmpty() ? 0 : spectra.ReverseBegin()->first + 1);

      // add the spectrum to the misc properties
      MOLECULE.StoreProperty( spectrum_name_base + util::Format()( count), value);

      // add properties

      // field strength
      storage::Pair< size_t, std::string> field_strength( count, s_UnreportedDescriptor);
      storage::Map< size_t, std::string> field_strength_props( GetPropertyMap( MOLECULE, s_FieldStrengthPropertyDescriptor));
      if( util::IsDefined( m_FieldStrength))
      {
        field_strength.Second() = util::Format()( m_FieldStrength);
      }
      BCL_Assert( field_strength_props.Insert( field_strength).second, "field strength already set: " + util::Format()( field_strength));
      MOLECULE.StoreProperty( s_FieldStrengthPropertyDescriptor, MapToLine( field_strength_props));

      // solvent
      storage::Map< size_t, std::string> solvent_props( GetPropertyMap( MOLECULE, s_SolventPropertyDescriptor));
      storage::Pair< size_t, std::string> solvent( count, s_UnreportedDescriptor);
      if( !m_Solvent.empty())
      {
        solvent.Second() = m_Solvent;
      }
      BCL_Assert( solvent_props.Insert( solvent).second, "solvent already set: " + util::Format()( solvent));
      MOLECULE.StoreProperty( s_SolventPropertyDescriptor, MapToLine( solvent_props));

      // assignment method
      storage::Map< size_t, std::string> assignment_method_props( GetPropertyMap( MOLECULE, s_AssignmentMethodDescriptor));
      storage::Pair< size_t, std::string> assignment( count, s_UnreportedDescriptor);
      if( !m_AssignmentMethod.empty())
      {
        assignment.Second() = m_AssignmentMethod;
      }
      BCL_Assert( assignment_method_props.Insert( assignment).second, "assignment already set: " + util::Format()( assignment));
      MOLECULE.StoreProperty( s_AssignmentMethodDescriptor, MapToLine( assignment_method_props));

      // temperature
      storage::Map< size_t, std::string> temperature_props( GetPropertyMap( MOLECULE, s_TemperaturePropertyDescriptor));
      storage::Pair< size_t, std::string> temperature( count, s_UnreportedDescriptor);
      if( util::IsDefined( m_Temperature))
      {
        temperature.Second() = util::Format()( m_Temperature);
      }
      BCL_Assert( temperature_props.Insert( temperature).second, "temperature already set: " + util::Format()( temperature));
      MOLECULE.StoreProperty( s_TemperaturePropertyDescriptor, MapToLine( temperature_props));
    }

    //! @brief get all signals involving given atom
    //! @param ATOM atom of interest
    //! @return SiPtrVector< const Signal> all signals involving atom of interest
    util::SiPtrVector< const Signal> Spectrum::DetermineInvolvedSignals( const chemistry::AtomConformationalInterface &ATOM) const
    {
      // all collected signals
      util::SiPtrVector< const Signal> signals;

      // iterate over all Signals in that Spectrum
      for
      (
        util::ShPtrVector< Signal>::const_iterator itr_signals( m_Signals.Begin()), itr_signals_end( m_Signals.End());
        itr_signals != itr_signals_end;
        ++itr_signals
      )
      {
        // check if current signal contains this atom
        if( ( *itr_signals)->ContainsAtom( ATOM))
        {
          signals.PushBack( util::SiPtr< const Signal>( **itr_signals));
        }
      }

      // return all signals for that atom
      return signals;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Spectrum from std::istream
    //! @param ISTREAM input stream that contains spectrum object
    std::istream &Spectrum::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_Type            , ISTREAM);
      io::Serialize::Read( m_Signals         , ISTREAM);
      io::Serialize::Read( m_Temperature     , ISTREAM);
      io::Serialize::Read( m_Solvent         , ISTREAM);
      io::Serialize::Read( m_FieldStrength   , ISTREAM);
      io::Serialize::Read( m_AssignmentMethod, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write Spectrum into std::ostream
    //! @param OSTREAM output stream in which the spectrum object is written
    //! @param INDENT indentation
    std::ostream &Spectrum::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_Type            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Signals         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Temperature     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Solvent         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FieldStrength   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AssignmentMethod, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert property line to property map
    //! @param CONFORMATION the conformation that should have the given property
    //! @param PROPERTY_LINE string of the property line "0:Unreported 1:CDCl3" create a map of size_t:string
    //! @return map that maps id to property
    storage::Map< size_t, std::string> Spectrum::GetPropertyMap
    (
      const chemistry::ConformationInterface &CONFORMATION,
      const std::string &PROPERTY_LINE
    )
    {
      const std::string &line( CONFORMATION.GetMDLProperty( PROPERTY_LINE));
      storage::Map< size_t, std::string> property_map;

      // split the string id separator
      const storage::Vector< std::string> split1
      (
        util::SplitString( line, std::string( 1, s_SpectraDescriptorSeparator))
      );

      if( split1.IsEmpty())
      {
        return property_map;
      }

      // get index of first element
      BCL_Assert
      (
        util::IsNumerical( split1.FirstElement()),
        "line needs to begin with an index, but is: " + line
      );
      size_t index( util::ConvertStringToNumericalValue< size_t>( split1.FirstElement()));
      for
      (
        storage::Vector< std::string>::const_iterator itr( split1.Begin() + 1), itr_last( split1.Last());
        itr <= itr_last;
        ++itr
      )
      {
        std::pair< size_t, std::string> key_value( index, "");

        // itr does not contain another index at the end
        if( itr == itr_last)
        {
          key_value.second = util::TrimString( *itr);
        }
        // itr contains the index for the next property at the end
        else
        {
          // find the property separator from the back
          const std::string::size_type split_pos( itr->rfind( s_PropertySeparator));
          BCL_Assert( split_pos != std::string::npos, "missing index for last property in line: " + line);

          // index is at the end
          const std::string next_index( itr->substr( split_pos, itr->length() - split_pos));
          BCL_Assert( util::IsNumerical( next_index), "non numerical index in line: " + line);
          index = util::ConvertStringToNumericalValue< size_t>( next_index);

          // actual property is before the index
          key_value.second = util::TrimString( itr->substr( 0, split_pos));
        }

        // check that it can be inserted without overlap of keys
        BCL_Assert( property_map.Insert( key_value).second, "properties with identical index: " + line);
      }

      // end
      return property_map;
    }

    //! @brief convert property map to property line "0:Unreported 1:CDCl3"
    //! @param PROPERTY_MAP that maps size_t to string
    //! @return string of the form "0:Unreported 1:CDCl3"
    std::string Spectrum::MapToLine( const storage::Map< size_t, std::string> &PROPERTY_MAP)
    {
      std::string property_line;

      // iterate over map
      for
      (
        storage::Map< size_t, std::string>::const_iterator itr( PROPERTY_MAP.Begin()), itr_end( PROPERTY_MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        property_line += util::Format()( itr->first);
        property_line += s_SpectraDescriptorSeparator;
        property_line += itr->second;
        property_line += s_PropertySeparator;
      }

      // end
      return property_line;
    }

    //! @brief read signals from and MDL Spectrum line
    //! @param MOLECULE the small molecule the spectrum is for
    //! @param SPECTRUM_LINE the spectrum line with all the signals from mdl file misc properties in NMRSHIFTDB
    //! @return a SpHptrVector of Signals
    util::ShPtrVector< Signal> Spectrum::ReadSignals
    (
      const chemistry::ConformationInterface &MOLECULE,
      const std::string &SPECTRUM_LINE
    )
    {
      // shared pointer vector of signal1Ds
      util::ShPtrVector< Signal> vector_signal;

      // split spectrum line into signals
      const storage::Vector< std::string> signal_strings( util::SplitString( SPECTRUM_LINE, std::string( 1, s_SignalSeparator)));

      // iterate over all signal string
      for
      (
        storage::Vector< std::string>::const_iterator itr( signal_strings.Begin()), itr_end( signal_strings.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->length() < 2)
        {
          continue;
        }

        vector_signal.PushBack
        (
          util::ShPtr< Signal>
          (
            new Signal( util::ShPtrVector< Signal1D>( 1, Signal1D::SignalFromString( *itr, MOLECULE.GetAtomsIterator())))
          )
        );
      }

      // end
      return vector_signal;
    }

    //! @brief generate the spectra
    //! @param MOLECULE the small molecule associated with the spectra
    //! @return a Map of spectra, where the key is the spectra number
    storage::Map< size_t, util::ShPtr< Spectrum> > Spectrum::GenerateSpectra
    (
      const chemistry::ConformationInterface &MOLECULE
    )
    {
      // spectra properties
      const storage::Map< size_t, std::string> field_strength_props(    GetPropertyMap( MOLECULE, Spectrum::s_FieldStrengthPropertyDescriptor));
      const storage::Map< size_t, std::string> solvent_props(           GetPropertyMap( MOLECULE, Spectrum::s_SolventPropertyDescriptor));
      const storage::Map< size_t, std::string> temperature_props(       GetPropertyMap( MOLECULE, Spectrum::s_TemperaturePropertyDescriptor));
      const storage::Map< size_t, std::string> assignment_method_props( GetPropertyMap( MOLECULE, Spectrum::s_AssignmentMethodDescriptor));

      // initialize return value
      storage::Map< size_t, util::ShPtr< Spectrum> > spectra;

      // iterate over all properties and find Spectrum line
      for
      (
        storage::Map< std::string, std::string>::const_iterator
          prop_itr( MOLECULE.GetStoredProperties().Begin()), prop_itr_end( MOLECULE.GetStoredProperties().End());
        prop_itr != prop_itr_end;
        ++prop_itr
      )
      {
        const storage::Vector< std::string> descriptors( util::SplitString( prop_itr->first));
        if( descriptors( 0) != s_SpectrumIdentifier)
        {
          continue;
        }

        //BCL_MessageStd( MOLECULE.GetMiscellaneousProperty( "Spectrum 13C 0"));
        const size_t current_number( util::ConvertStringToNumericalValue< size_t>( descriptors( 2)));

        // check that spectra id is unique
        if( spectra.Find( current_number) != spectra.End())
        {
          BCL_MessageCrt( "spectra with duplicate id: " + descriptors( 2));
          continue;
        }

        // read out spectrum line in MDL block if available
        const util::ShPtrVector< Signal> signals
        (
          Spectrum::ReadSignals( MOLECULE, prop_itr->second)
        );

        // spectrum
        util::ShPtr< Spectrum> current_spectrum
        (
          new Spectrum( Spectrum::SpecTypeEnum( descriptors( 1)), signals)
        );

        // set temp and field strength, solvent and assignment method
        const storage::Map< size_t, std::string>::const_iterator temp_itr( temperature_props.Find(       current_number));
        const storage::Map< size_t, std::string>::const_iterator fiel_itr( field_strength_props.Find(    current_number));
        const storage::Map< size_t, std::string>::const_iterator solv_itr( solvent_props.Find(           current_number));
        const storage::Map< size_t, std::string>::const_iterator assi_itr( assignment_method_props.Find( current_number));

        // temperature
        if( temp_itr == temperature_props.End() || !util::IsNumerical( temp_itr->second))
        {
          // no temperature given
          current_spectrum->SetTemperature( 298.15);
        }
        else
        {
          current_spectrum->SetTemperature( util::ConvertStringToNumericalValue< double>( temp_itr->second));
        }

        // field strength
        if( fiel_itr == field_strength_props.End() || !util::IsNumerical( fiel_itr->second))
        {
          // no field strength given
          current_spectrum->SetFieldStrength( util::GetUndefined< double>());
        }
        else
        {
          current_spectrum->SetFieldStrength( util::ConvertStringToNumericalValue< double>( fiel_itr->second));
        }

        // solvent
        if( solv_itr == solvent_props.End() || solv_itr->second.empty())
        {
          current_spectrum->SetSolvent( Spectrum::s_UnreportedDescriptor);
        }
        else
        {
          current_spectrum->SetSolvent( solv_itr->second);
        }

        // assignment method
        if( assi_itr == assignment_method_props.End() || assi_itr->second.empty())
        {
          current_spectrum->SetAssignmentMethod( Spectrum::s_UnreportedDescriptor);
        }
        else
        {
          current_spectrum->SetAssignmentMethod( assi_itr->second);
        }

        // insert spectrum into map
        spectra[ current_number] = current_spectrum;
      }

      // end
      return spectra;
    }

  } // namespace nmr
} // namespace bcl
