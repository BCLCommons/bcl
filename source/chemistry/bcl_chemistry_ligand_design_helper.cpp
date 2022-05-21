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
#include "chemistry/bcl_chemistry_ligand_design_helper.h"

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "model/bcl_model.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_mdl_handler.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LigandDesignHelperImpl
    //! @brief Functionality needed by the (java-based) BCL::LigandDesign tool
    //!
    //! @see @link example_chemistry_ligand_design_helper.cpp @endlink
    //! @author mendenjl
    //! @date Mar 02, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LigandDesignHelperImpl :
      public util::SerializableInterface
    {

    private:

      SampleConformations m_ConformationSampler; //!< Generates a conformation for a designed molecule

      //! All the properties the user can choose from
      //! Stored as a comma-separated string for simplified java inter-exchange
      std::string         m_AllPropertiesCsv;

      //! All the predictions the user can choose from
      //! Stored as a comma-separated string for simplified java inter-exchange
      std::string         m_AllPredictionsCsv;

      //! Map from property name to implementation
      storage::Map< std::string, descriptor::CheminfoProperty> m_PropertyMap;

      //! Map from property name to implementation
      storage::Map< std::string, descriptor::CheminfoProperty> m_PredictionMap;

      //! directory containing the jar file that uses this class
      std::string m_JarFileDirectory;

      //! Mutex to prevent multiple java threads from using the descriptors at the same time
      sched::Mutex m_Mutex;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LigandDesignHelperImpl();

      //! @brief constructor from jar file directory path
      explicit LigandDesignHelperImpl( const std::string &JAR_FILE_DIRECTORY);

      //! virtual copy constructor
      LigandDesignHelperImpl *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief get the comma separated properties list
      //! @return the comma-seperated properties list
      const std::string &GetAvailableProperties() const;

      //! @brief get the comma separated predictions list
      //! @note prediction is defined here as anything that requires a conformation
      //! @return the comma-seperated predictions list
      const std::string &GetAvailablePredictions() const;

      //! @brief set the directory of the jar file. This is used to locate the bcl-related files
      //! @param JAR_FILE_DIRECTORY the directory containing the jar file; should also contain a bcl folder
      void SetJarDirectory( const std::string &JAR_FILE_DIR);

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate conformation-independent properties.
      //! @param MOLECULE molecule of interest, encoded as a string in SDF format
      //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
      //! @return DESCRIPTORS given as values in a comma separated list.
      //!         Elements in vector-valued descriptors separated by spaces
      //!         Note that if the atom types are undefined, all values will be "Unknown"
      std::string CalculateProperties
      (
        const std::string &MOLECULE,
        const std::string &PROPERTIES
      );

      //! @brief processes the molecule; producing a 3d-conformation,
      //! @param MOLECULE molecule of interest, encoded as a string in SDF format
      //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
      //! @return molecule in sdf formation with 3D conformation and DESCRIPTORS as misc properties
      std::string Process
      (
        const std::string &MOLECULE,
        const std::string &PROPERTIES
      );

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LigandDesignHelperImpl::LigandDesignHelperImpl() :
      m_ConformationSampler( RotamerLibraryFile(), "bcl::chemistry::ConformationComparisonInterface", 0.0, 5, 10, false),
      m_AllPropertiesCsv(),
      m_PropertyMap(),
      m_PredictionMap(),
      m_JarFileDirectory()
    {
    }

    //! @brief constructor from jar file directory path
    LigandDesignHelperImpl::LigandDesignHelperImpl( const std::string &JAR_FILE_DIRECTORY) :
      m_ConformationSampler( RotamerLibraryFile(), "bcl::chemistry::ConformationComparisonInterface", 0.0, 5, 10, false),
      m_AllPropertiesCsv(),
      m_PropertyMap(),
      m_PredictionMap(),
      m_JarFileDirectory()
    {
      SetJarDirectory( JAR_FILE_DIRECTORY);
    }

    //! virtual copy constructor
    LigandDesignHelperImpl *LigandDesignHelperImpl::Clone() const
    {
      return new LigandDesignHelperImpl( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &LigandDesignHelperImpl::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &LigandDesignHelperImpl::GetAlias() const
    {
      static const std::string s_name( "LigandDesignHelperImpl");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief get the comma separated properties list
    //! @return the comma-seperated properties list
    const std::string &LigandDesignHelperImpl::GetAvailableProperties() const
    {
      return m_AllPropertiesCsv;
    }

    //! @brief get the comma separated predictions list
    //! @note prediction is defined here as anything that requires a conformation
    //! @return the comma-seperated predictions list
    const std::string &LigandDesignHelperImpl::GetAvailablePredictions() const
    {
      return m_AllPredictionsCsv;
    }

    //! @brief set the directory of the jar file. This is used to locate the bcl-related files
    //! @param JAR_FILE_DIRECTORY the directory containing the jar file; should also contain a bcl folder
    void LigandDesignHelperImpl::SetJarDirectory( const std::string &JAR_FILE_DIR)
    {
      m_Mutex.Lock();
      if( m_JarFileDirectory == JAR_FILE_DIR)
      {
        m_Mutex.Unlock();
        return;
      }
      m_JarFileDirectory = JAR_FILE_DIR;

      // make any missing directories
      io::Directory::MkDir( JAR_FILE_DIR + "/bcl");
      io::Directory::MkDir( JAR_FILE_DIR + "/bcl/qsar");

      // set the model directory
      model::Model::GetModelPathFlag()->ReadFromList
      (
        storage::Vector< std::string>::Create( JAR_FILE_DIR + "/bcl/qsar"), util::GetLogger()
      );
      app::GetApps().GetExecutablePath() = JAR_FILE_DIR + "/bcl";

      // write an example configuration if none exists
      const std::string config_filename( JAR_FILE_DIR + "/bcl/config.txt");

      if( io::DirectoryEntry( config_filename).DoesExist())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, config_filename);
        BCL_Assert
        (
          this->TryRead( util::ObjectDataLabel( input), util::GetLogger()),
          "Couldn't read configuration file at " + config_filename
        );
        io::File::CloseClearFStream( input);
      }
      else
      {
        // write an example config file
        m_PropertyMap[ "Weight"] = descriptor::GetCheminfoProperties().calc_Mass;
        m_PropertyMap[ "NRotBond"] = descriptor::GetCheminfoProperties().calc_NRotBond;
        io::OFStream output;
        io::File::MustOpenOFStream( output, config_filename);
        output << this->GetLabel().ToStringDefaultWidth();
        io::File::CloseClearFStream( output);
      }

      // load in the property map and sample conformations
      io::IFStream input;
      io::File::MustOpenIFStream( input, config_filename);
      this->AssertRead( util::ObjectDataLabel( input));
      io::File::CloseClearFStream( input);

      // update m_AllPropertiesCsv based on property map
      storage::Set< std::string> properties_set( m_PropertyMap.GetKeys());
      m_AllPropertiesCsv = util::Join( ",", storage::Vector< std::string>( properties_set.Begin(), properties_set.End()));
      storage::Set< std::string> predictions_set( m_PredictionMap.GetKeys());
      m_AllPredictionsCsv = util::Join( ",", storage::Vector< std::string>( predictions_set.Begin(), predictions_set.End()));
      m_Mutex.Unlock();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate conformation-independent properties.
    //! @param MOLECULE molecule of interest, encoded as a string in SDF format
    //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
    //! @return DESCRIPTORS given as values in a comma separated list.
    //!         Elements in vector-valued descriptors separated by spaces
    //!         Note that if the atom types are undefined, returns empty string
    std::string LigandDesignHelperImpl::CalculateProperties
    (
      const std::string &MOLECULE,
      const std::string &PROPERTIES
    )
    {
      m_Mutex.Lock();
      std::stringstream str( MOLECULE);
      sdf::MdlHandler handler( str);
      FragmentComplete frag( sdf::FragmentFactory::MakeFragment( handler, sdf::e_Saturate));
      if( frag.HasNonGasteigerAtomTypes())
      {
        m_Mutex.Unlock();
        BCL_MessageStd( "Non-gasteiger atom types: " + frag.GetAtomTypesString());
        return "";
      }

      storage::Vector< std::string> properties( util::SplitString( PROPERTIES, ","));
      std::ostringstream output;
      for
      (
        storage::Vector< std::string>::const_iterator itr( properties.Begin()), itr_end( properties.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Map< std::string, descriptor::CheminfoProperty>::iterator itr_map_prop
        (
          m_PropertyMap.Find( *itr)
        );
        BCL_Assert
        (
          itr_map_prop != m_PropertyMap.End(),
          "Given property or prediction: " + *itr + " is not available; available properties are: " + m_AllPropertiesCsv
        );
        linal::Vector< float> prop_values( itr_map_prop->second->SumOverObject( frag));
        if( itr != properties.Begin())
        {
          output << ", ";
        }
        for
        (
          linal::Vector< float>::const_iterator itr_prop( prop_values.Begin()), itr_prop_end( prop_values.End());
          itr_prop != itr_prop_end;
          ++itr_prop
        )
        {
          output << *itr_prop;
          if( itr_prop + 1 != itr_prop_end)
          {
            output << ' ';
          }
        }
      }
      m_Mutex.Unlock();
      return output.str();
    }

    //! @brief processes the molecule; producing a 3d-conformation,
    //! @param MOLECULE molecule of interest, encoded as a string in SDF format
    //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
    //! @return molecule in sdf formation with 3D conformation and DESCRIPTORS as misc properties
    std::string LigandDesignHelperImpl::Process
    (
      const std::string &MOLECULE,
      const std::string &PREDICTIONS
    )
    {
      m_Mutex.Lock();
      std::stringstream str( MOLECULE);
      sdf::MdlHandler handler( str);
      FragmentComplete frag( sdf::FragmentFactory::MakeFragment( handler, sdf::e_Saturate));
      if( frag.HasNonGasteigerAtomTypes())
      {
        m_Mutex.Unlock();
        BCL_MessageStd( "Non-gasteiger atom types: " + frag.GetAtomTypesString());
        return std::string();
      }
      random::GetGlobalRandom().SetGlobalSeedFromCommandlineFlag();
      FragmentEnsemble ensemble( m_ConformationSampler( frag).First());
      if( ensemble.GetSize() == size_t( 0))
      {
        m_Mutex.Unlock();
        BCL_MessageStd( "Conformation could not be generated");
        return std::string();
      }
      else
      {
        frag = ensemble.GetMolecules().FirstElement();
      }

      storage::Vector< std::string> properties( util::SplitString( PREDICTIONS, ","));
      for
      (
        storage::Vector< std::string>::const_iterator itr( properties.Begin()), itr_end( properties.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Map< std::string, descriptor::CheminfoProperty>::iterator itr_map_pred
        (
          m_PredictionMap.Find( *itr)
        );
        if( itr_map_pred != m_PredictionMap.End())
        {
          frag.StoreProperty( *itr, itr_map_pred->second->SumOverObject( frag));
          continue;
        }
        storage::Map< std::string, descriptor::CheminfoProperty>::iterator itr_map_prop
        (
          m_PropertyMap.Find( *itr)
        );
        BCL_Assert
        (
          itr_map_prop != m_PropertyMap.End(),
          "Given property or prediction: " + *itr + " is not available; available properties are: " + m_AllPropertiesCsv
          + " available predictions are: " + m_AllPredictionsCsv
        );
        frag.StoreProperty( *itr, itr_map_prop->second->SumOverObject( frag));
      }
      std::stringstream outs;
      frag.WriteMDL( outs);
      m_Mutex.Unlock();
      return outs.str();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LigandDesignHelperImpl::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.AddInitializer
      (
        "properties",
        "map of desired property name to actual descriptor initializer",
        io::Serialization::GetAgent( &m_PropertyMap)
      );
      serializer.AddInitializer
      (
        "predictions",
        "map of desired prediction name to actual descriptor initializer",
        io::Serialization::GetAgent( &m_PredictionMap)
      );
      serializer.AddInitializer
      (
        "conformation creator",
        "parameters for sampling the best conformation",
        io::Serialization::GetAgent( &m_ConformationSampler)
      );

      return serializer;
    }

    //! @brief singleton retriever class for ligand design helper implementation
    LigandDesignHelperImpl &GetLigandDesignHelperImpl()
    {
      static LigandDesignHelperImpl s_impl;
      return s_impl;
    }

    //! @brief get the comma separated properties list
    //! @return the comma-seperated properties list
    const std::string &LigandDesignHelper::GetAvailableProperties()
    {
      return GetLigandDesignHelperImpl().GetAvailableProperties();
    }

    //! @brief get the comma separated properties list
    //! @return the comma-seperated properties list
    const std::string &LigandDesignHelper::GetAvailablePredictions()
    {
      return GetLigandDesignHelperImpl().GetAvailablePredictions();
    }

    //! @brief set the directory of the jar file. This is used to locate the bcl-related files
    //! @param JAR_FILE_DIRECTORY the directory containing the jar file; should also contain a bcl folder
    void LigandDesignHelper::SetJarDirectory( const std::string &JAR_FILE_DIR)
    {
      return GetLigandDesignHelperImpl().SetJarDirectory( JAR_FILE_DIR);
    }

    //! @brief calculate conformation-independent properties.
    //! @param MOLECULE molecule of interest, encoded as a string in SDF format
    //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
    //! @return DESCRIPTORS given as values in a comma separated list.
    //!         Elements in vector-valued descriptors separated by spaces
    //!         Note that if the atom types are undefined, all values will be "Unknown"
    std::string LigandDesignHelper::CalculateProperties
    (
      const std::string &MOLECULE,
      const std::string &PROPERTIES
    )
    {
      return GetLigandDesignHelperImpl().CalculateProperties( MOLECULE, PROPERTIES);
    }

    //! @brief processes the molecule; producing a 3d-conformation,
    //! @param MOLECULE molecule of interest, encoded as a string in SDF format
    //! @param DESCRIPTORS descriptors (comma seperated) to calculate for the molecule
    //! @return molecule in sdf formation with 3D conformation and DESCRIPTORS as misc properties
    std::string LigandDesignHelper::ProcessMolecule
    (
      const std::string &MOLECULE,
      const std::string &DESCRIPTORS
    )
    {
      return GetLigandDesignHelperImpl().Process( MOLECULE, DESCRIPTORS);
    }

  } // namespace chemistry
} // namespace bcl
