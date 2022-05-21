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
#include "bcl_app_molecule_properties.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_string_property_interface.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_default.h"
#include "command/bcl_command_parameter_check_or.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_base.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_string_replacement.h"

namespace bcl
{
  namespace app
  {
  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    MoleculeProperties::MoleculeProperties() :
      m_HistogramNumericPropertyFlag
      (
        new command::FlagDynamic
        (
          "numeric_histogram",
          "properties containing numeric values, min value, bin size, # bins",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "property",
              "atom or small molecule property",
              command::ParameterCheckSerializable( descriptor::CheminfoProperty())
            ),
            command::Parameter
            (
              "min_value",
              "min/left boundary of the histogram",
              command::ParameterCheckRanged< double>()
            ),
            command::Parameter
            (
              "bin_size",
              "size of each bin in the histogram",
              command::ParameterCheckRanged< double>( 0, std::numeric_limits< double>::max())
            ),
            command::Parameter
            (
              "number_bins",
              "number of bins in the histogram",
              command::ParameterCheckRanged< size_t>()
            )
          )
        )
      ),
      m_HistogramStringPropertyFlag
      (
        new command::FlagDynamic
        (
          "string_histogram",
          "properties containing space-delimited strings",
          command::Parameter
          (
            "property",
            "property",
            command::ParameterCheckSerializable( util::Implementation< chemistry::StringPropertyInterface>())
          )
        )
      ),
      m_HistogramOutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_histogram",
          "filename to write out histograms",
          command::Parameter
          (
            "output_histogram",
            "filename to write out histograms",
            "molecule.histogram.txt"
          )
        )
      ),
      m_AveMinMaxStdPropertyFlag
      (
        new command::FlagDynamic
        (
          "statistics",
          "properties on which to take statistics (output file is controlled by -output_histogram)",
          command::Parameter
          (
            "statistics",
            "properties on which to take statistics (output file is controlled by -output_histogram)",
            command::ParameterCheckSerializable( util::Implementation< chemistry::StringPropertyInterface>())
          )
        )
      ),
      m_TablePropertiesFlag
      (
        new command::FlagDynamic
        (
          "tabulate",
          "properties to put into a csv file; file will always contain molecule index in the first column",
          command::Parameter
          (
            "property",
            "property",
            command::ParameterCheckSerializable( util::Implementation< chemistry::StringPropertyInterface>())
          )
        )
      ),
      m_TableOutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_table",
          "filename to write out csv of properties given in -table",
          command::Parameter
          (
            "output_table",
            "filename to write out csv of properties given in -table",
            "properties.csv"
          )
        )
      ),
      m_RemovePropertiesFlag
      (
        new command::FlagDynamic
        (
          "remove",
          "miscellaneous properties to remove from all molecules",
          command::Parameter
          (
            "property",
            "miscellaneous property to remove from all molecules"
          )
        )
      ),
      m_RemoveAllPropertiesFlag
      (
        new command::FlagStatic
        (
          "remove_all",
          "removes all properties from all molecules"
        )
      ),
      m_RenamePropertiesFlag
      (
        new command::FlagDynamic
        (
          "rename",
          "for each pair of properties listed, change the name from the first property to the second",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "initial_property_name",
              "property name as it is currently in the molecule"
            ),
            command::Parameter
            (
              "desired_property_name",
              "property name to rename the initial property name to"
            )
          )
        )
      ),
      m_AddPropertiesFlag
      (
        new command::FlagDynamic
        (
          "add",
          "descriptors / properties to add to all molecules",
          command::Parameter
          (
            "property",
            "descriptors / properties to add to all molecules",
            command::ParameterCheckOr
            (
              command::ParameterCheckSerializable( util::Implementation< chemistry::StringPropertyInterface>()),
              command::ParameterCheckAllowed( storage::Vector< std::string>( 1, "Index")),
              command::ParameterCheckSerializable( descriptor::CheminfoProperty())
            )
          )
        )
      ),
      m_AddPropertyStringsFlag
      (
        new command::FlagDynamic
        (
          "add_strings",
          "strings to add as descriptors / properties to all molecules",
          storage::Vector< command::Parameter>::Create
          (
            command::Parameter
            (
              "property",
              "name of descriptor / property to add"
            ),
            command::Parameter
            (
              "string",
              "the string assigned to the property / descriptor name for all molecules"
            )
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "sdf filename for where to write out molecules",
          command::Parameter
          (
            "filename",
            "sdf filename for where to write out molecules",
            command::ParameterCheckDefault(),
            ""
          )
        )
      )
    {
    }

    //! @brief copy constructor; only copy flags for applications
    MoleculeProperties::MoleculeProperties( const MoleculeProperties &PARENT) :
      m_HistogramNumericPropertyFlag( PARENT.m_HistogramNumericPropertyFlag),
      m_HistogramStringPropertyFlag( PARENT.m_HistogramStringPropertyFlag),
      m_HistogramOutputFilenameFlag( PARENT.m_HistogramOutputFilenameFlag),
      m_AveMinMaxStdPropertyFlag( PARENT.m_AveMinMaxStdPropertyFlag),
      m_TablePropertiesFlag( PARENT.m_TablePropertiesFlag),
      m_TableOutputFilenameFlag( PARENT.m_TableOutputFilenameFlag),
      m_RemovePropertiesFlag( PARENT.m_RemovePropertiesFlag),
      m_RemoveAllPropertiesFlag( PARENT.m_RemoveAllPropertiesFlag),
      m_RenamePropertiesFlag( PARENT.m_RenamePropertiesFlag),
      m_AddPropertiesFlag( PARENT.m_AddPropertiesFlag),
      m_AddPropertyStringsFlag( PARENT.m_AddPropertyStringsFlag),
      m_OutputFilenameFlag( PARENT.m_OutputFilenameFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeProperties *MoleculeProperties::Clone() const
    {
      return new MoleculeProperties( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeProperties::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeProperties::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // get the data source
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // numeric property to take basic histogram of
      sp_cmd->AddFlag( m_HistogramNumericPropertyFlag);

      // string property to take histogram of
      sp_cmd->AddFlag( m_HistogramStringPropertyFlag);

      // flag for output filename for histograms
      sp_cmd->AddFlag( m_HistogramOutputFilenameFlag);

      // flag for properties on which to take statistics
      sp_cmd->AddFlag( m_AveMinMaxStdPropertyFlag);

      // properties from which to make a table
      sp_cmd->AddFlag( m_TablePropertiesFlag);

      // filename to write out property table to
      sp_cmd->AddFlag( m_TableOutputFilenameFlag);

      // any properties to remove
      sp_cmd->AddFlag( m_RemovePropertiesFlag);

      // remove all properties to remove
      sp_cmd->AddFlag( m_RemoveAllPropertiesFlag);

      // any properties that should be added
      sp_cmd->AddFlag( m_AddPropertiesFlag);

      // any strings that should be added as properties
      sp_cmd->AddFlag( m_AddPropertyStringsFlag);

      // any properties to rename
      sp_cmd->AddFlag( m_RenamePropertiesFlag);

      // output for molecules, if properties were changed
      sp_cmd->AddFlag( m_OutputFilenameFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeProperties::GetDescription() const
    {
      return "Works with string and numeric descriptions of molecules";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeProperties::GetReadMe() const
    {
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of the BCL application molecule:Properties, its terms of use, "
        "appropriate citation, installation procedures, molecule:Properties execution, and technical support\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. What is molecule:Properties?\n"
        "molecule:Properties is a C++ based application, created by Vanderbilt University's Meiler Laboratory, "
        "which is part a larger library of applications called BCL::Commons. molecule:Properties works with "
        "numeric and string descriptions (called properties) of molecules. Properties can be added to an SDF file, "
        "renamed, or removed.  Statistics and histograms can also be collected for properties.  Lastly, property values "
        " can be tabulated for analysis with external programs."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING molecule:Properties.\n"
        "When using descriptor:Analyze in a publication, please cite the following publications describing the "
        "application's development:\n"
        "\n"
        "Butkiewicz M, Lowe EW, Mueller R, Mendenhall JL, Teixeira PL, Weaver CD, Meiler J. "
        "Benchmarking Ligand-Based Virtual High-Throughput Screening with the PubChem Database "
        "Molecules, 18, (1), 735-756. ; 2013\n"
        "Link:  www.http://meilerlab.org/index.php/publications/show/2013\nJournal link: http://www.mdpi.com/1420-3049/18/1/735\n"
        "\n"
        "If the MolecularAsymmetry descriptor is used, the following publication must also be cited: "
        "Sliwoski, Gregory, et al. \"BCL:: EMAS—Enantioselective Molecular Asymmetry Descriptor for 3D-QSAR.\" "
        "Molecules 17.8 (2012): 9971-9989.\n"
        "Link:  www.http://meilerlab.org/index.php/publications/show/2012\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING molecule:Properties.\n"
        "molecule:Properties can be run in either or both of two primary modes: property manipulation or analysis\n"
        "Here are some ways to learn what properties are available and how to work with them:\n"
        "bcl.exe molecule:Properties -help\n  shows lists of all properties for each flag, individually\n"
        "bcl.exe molecule:Properties -statistics help shows lists of all properties that can be used with the statistics command\n"
        "VI. A. Using conventional, non-parametric properties\n"
        "Most properties, such as Atom_SigmaCharge, NRings, NAtoms, etc. are not parameterized or user customizable, and "
        "so can be used without further specification, e.g.\nbcl.exe -input_filenames ./molecules.sdf -output_filenames"
        "./molecules.wprops.sdf -add NAtoms NAromaticRings HBondAcceptor\n"
        "VI. B. Creating novel properties\n"
        "The properties framework is highly flexible and allows user specification of novel meta-properties. Properties"
        " that can be used to create other novel property combinations using a function-call like syntax known as OGDL ("
        "Ordered Graph Data Language, see http://ogdl.org/ for details).  However, detailed knowledge of OGDL is "
        "unnecessary for most use cases because help is always available from the bcl.  For example, running\n"
        "bcl.exe molecule:Properties -statistics help will show all numeric properties, among them is the following line:\n"
        "* RDF : computes the radial distribution function using a given atom property; \"RDF(help)\" shows internal options\n"
        "which indicates that one can receive help by typing the following:\nbcl.exe molecule:Properties -statistics 'RDF(help)'\n"
        "which yields:\n"
        "   computes the radial distribution function using a given atom property\n"
        "   Default label : RDF(property=Atom_Identity,step size=0.1,temperature=100,steps=128,normalized=0)\n"
        "   Parameters:\n"
        "   <property> property over which to calculate the radial distribution function,   1 for every atom\n"
        "   <step size> size of each step in angstroms, default: \"0.25\", Any decimal (floating-point) value >= 0.01 <= 100\n"
        "   <temperature> increasing temperature spreads autocorrelation across more distant bins, default: \"100\", Any decimal (floating-point) value >= 0 <= 1000\n"
        "   <steps> # of steps/bins (each of size = step size) used in the radial distribution function, default: \"48\", Any non-negative integer >= 1 <= 1000000\n"
        "   <normalized> whether to normalize the RDF by RDF of atom identity (1 for each atom), default: \"0\", Any non-negative integer\n"
        "all the parameters can also accept help, so to find out what you could use instead of Atom_Identity, execute:\n"
        "bcl.exe molecule:Properties -statistics 'RDF(property=help), which would then list all the atom properties that could be used\n"
        "Properties can also be multiplied, added, subtracted, concatenated, summed, logged, etc. using the appropriate meta-property\n"
        "Still other properties can load data from external files (Mapped), compute machine learning model predictions (PredictedActivity), and even retrieval of values from a database\n"
        "VI. C. SDF file property manipulation (Addition, Removal, Renaming)\n"
        "Property manipulation consists of up to three operations, which occur in the following order:\n"
        "1. adding properties to the sdf file specified with the -add flag\n"
        "2. renaming properties in the sdf file specified using the -rename flag\n"
        "3. removing properties in the sdf file specified using the -remove flag\n"
        "This ordering allows one to perform simultaneous addition, following by renaming of properties, e.g.\n"
        "bcl.exe -input_filenames molecules.sdf -output molecules_with_haccdon.sdf -add 'Add(HDon,HAcc)' -rename 'Add(HDon,HAcc)' HDonAcc\n"
        "which would effectively add a miscellaneous property HDonAcc to the sdf file, which is the sum of hydrogen bond donors and acceptors"
        "\n"
        "VI. D. Property analysis\n"
        "Properties can be analyzed via numeric or string histograms, statistics, and tabulation for external analysis\n"
        "see bcl.exe molecule:Properties -help for details.  Here are some example command lines: "
        "bcl.exe -input_filenames molecules.sdf -string_histogram AtomTypes BondTypes 'Numeric(NAtoms)' -output_histogram\n"
        "bcl.exe -input_filenames molecules.sdf -tabulate AtomTypes BondTypes 'Numeric(NAtoms)' -output_table table.txt\n"
        "bcl.exe -input_filenames molecules.sdf more_mols.sdf.bz2 -input_start 5 -input_max 50 -remove_h -numeric_histogram '2DA(property=Atom_SigmaCharge,steps=10)' -10 1 20 -output_histogram hist.txt\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF molecule:Properties.\n"
        "descriptor:Analyze is under ongoing further development. For current research please refer to "
        "www.meilerlab.org and navigate to research\n"
        + DefaultSectionSeparator()
      );

      return s_readme;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &MoleculeProperties::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL molecule:Properties works with numeric and string descriptions (called properties) of molecules.\n\n"
        "Features of BCL::MoleculeProperties\n"
        "<ul>"
        "  <li>Works with conventional properties like weight and total charge, as well as novel BCL descriptors like "
        "      smoothed, sign sensitive, 3D autocorrelations !molecule_properties.txt!</li>"
        "  <li>Gasteiger atom types, aromaticity and ring-aware bond-types, and chirality information can also be computed</li>"
        "  <li>Properties can be"
        "  <ul>\n"
        "    <li>added, removed, or renamed to an ensemble of molecules SDF file</li>\n"
        "    <li>Summarized with statistics (Mean/SD/Min/Max) over an ensemble</li>\n"
        "    <li>Summarized by histograms over an ensemble</li>\n"
        "    <li>Tabulated for analysis with external programs.</li>\n"
        "  </ul></li>\n\n"
        "  <li>Compressed molecule files (bz2, gzip) are supported</li>\n"
        "</ul>\n\n"
        "Example command line\n"
        "bcl.exe molecule:Properties -input_filenames 1_3_pentadiene_E.sdf -add AtomTypes BondTypes LogP Weight"
        " \"2DA(property=Atom_SigmaCharge)\" -output 1_3_pentadiene_E_with_props.sdf\n\n"
        "Output\n<code>\n"
        "Written by BCL::WriteToMDL,62204                                                           "
        "  -OEChem-05141111443D                                                                     "
        "                                                                                           "
        "  5  4  0  0  0  0  0  0  0  0999 V2000                                                    "
        "   -1.1769   -0.4059   -0.0001 C   0  0  0  0  0  3  0  0  0  0  0  0                      "
        "   -2.5154    0.2568    0.0001 C   0  0  0  0  0  4  0  0  0  0  0  0                      "
        "   -0.0151    0.2640   -0.0001 C   0  0  0  0  0  3  0  0  0  0  0  0                      "
        "    1.2729   -0.3894    0.0001 C   0  0  0  0  0  3  0  0  0  0  0  0                      "
        "    2.4344    0.2744    0.0000 C   0  0  0  0  0  3  0  0  0  0  0  0                      "
        "  1  2  1  0  0  0  0                                                                      "
        "  1  3  2  0  0  0  0                                                                      "
        "  3  4  1  0  0  0  0                                                                      "
        "  4  5  2  0  0  0  0                                                                      "
        "M  BCL ATM C_TrTrTrPi C_TeTeTeTe C_TrTrTrPi C_TrTrTrPi C_TrTrTrPi                          "
        "M  END                                                                                     "
        "> <2DA(property=Atom_SigmaCharge)>                                                         "
        "0.00156101 -0.000645298 -0.000125459 -9.44797e-06 -2.99351e-07 0 0 0 0 0 0                 "
        "                                                                                           "
        "> <AtomTypes>                                                                              "
        "C_TrTrTrPi C_TeTeTeTe C_TrTrTrPi C_TrTrTrPi C_TrTrTrPi                                     "
        "                                                                                           "
        "> <BondTypes>                                                                              "
        "NonConjugatedSingleBond ConjugatedDoubleBond_E ConjugatedSingleBond ConjugatedDoubleBond   "
        "                                                                                           "
        "> <LogP>                                                                                   "
        "1.51611                                                                                    "
        "                                                                                           "
        "> <Weight>                                                                                 "
        "60.05                                                                                      "
        "                                                                                           "
        "$$$$                                                                                       "
        "</code>"
      );

      return s_web_text;
    }

    //! @brief check that all the parameter choices were valid
    //! @return true if all the parameter choices were valid
    bool MoleculeProperties::CheckParametersAreAcceptable() const
    {
      if( m_TableOutputFilenameFlag->GetFlag() != m_TablePropertiesFlag->GetFlag())
      {
        BCL_MessageCrt
        (
          "Table properties and table output filename must be given together or not at all!"
        );
        return false;
      }
      if
      (
        m_HistogramOutputFilenameFlag->GetFlag() !=
        (
          m_HistogramNumericPropertyFlag->GetFlag()
          || m_HistogramStringPropertyFlag->GetFlag()
          || m_AveMinMaxStdPropertyFlag->GetFlag()
        )
      )
      {
        BCL_MessageCrt
        (
          "Histograms/Statistics and histogram output filename must be given together or not at all!"
        );
        return false;
      }
      if( !m_OutputFilenameFlag->GetFlag() && m_RemovePropertiesFlag->GetFlag())
      {
        BCL_MessageCrt
        (
          "Removing property has no effect unless given with an output filename!"
        );
      }

      return true;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeProperties::Main() const
    {
      if( !CheckParametersAreAcceptable())
      {
        return 0;
      }

      // Initialize the histograms and properties from the command line flags
      Initialize();

      // for all the files that were given on the command line
      for( chemistry::FragmentFeed itr_fragments; itr_fragments.NotAtEnd(); ++itr_fragments)
      {
        chemistry::FragmentComplete fragment( *itr_fragments);
        AddData( fragment, itr_fragments.GetPosition());
      }

      if( m_TableOutputFilenameFlag->GetFlag())
      {
        // close the output file
        io::File::CloseClearFStream( m_TableOutputFile);
      }
      if( m_OutputFilenameFlag->GetFlag())
      {
        // close the molecule output file
        io::File::CloseClearFStream( m_MoleculeOutputFile);
      }

      // open an output file for histograms
      if( m_HistogramOutputFilenameFlag->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          m_TableOutputFile,
          m_HistogramOutputFilenameFlag->GetFirstParameter()->GetValue()
        );
        // make histograms for each string property
        for
        (
          size_t property_number( 0), number_properties( m_StringHistograms.GetSize());
          property_number < number_properties;
          ++property_number
        )
        {
          // write the counts out to a string stream, which will then be written out to a message
          // this is to avoid having the =std=bcl::app> prefix on each line
          std::stringstream stream;
          stream << '\n'
                 << m_HistogramStringPropertyFlag->GetParameterList()( property_number)->GetValue() << " histogram\n"
                 << "Index\t\tCategory\t\t\tCounts\n";

          size_t index( 0);
          for
          (
            storage::Map< std::string, size_t>::const_iterator
              itr_prop( m_StringHistograms( property_number).Begin()),
              itr_prop_end( m_StringHistograms( property_number).End());
            itr_prop != itr_prop_end;
            ++itr_prop, ++index
          )
          {
            stream << index << "\t\t" << itr_prop->first << "\t\t" << itr_prop->second << '\n';
          }
          m_TableOutputFile << stream.str();
        }

        // make histograms for each numeric property
        for
        (
          size_t property_number( 0), number_properties( m_NumericHistograms.GetSize());
          property_number < number_properties;
          ++property_number
        )
        {
          // write the name of the property out (every 4th value in the command line flag)
          m_TableOutputFile
              << "\n"
              << m_HistogramNumericPropertyFlag->GetParameterList()( property_number * 4)->GetValue()
              << " histogram";

          // write out the histogram
          m_NumericHistograms( property_number).WriteVertically( m_TableOutputFile);
        }

        // write statistics for each property desired
        if( m_MeanStds.GetSize())
        {
          m_TableOutputFile << "\nProperty\tAve\tStd\tMin\tMax\tMax-Min\t(Max-Min)/Std";
        }
        for
        (
          size_t property_number( 0), number_properties( m_MeanStds.GetSize());
          property_number < number_properties;
          ++property_number
        )
        {
          const float stdev( m_MeanStds( property_number).GetStandardDeviation());
          const float range( m_MinMaxs( property_number).GetMax() - m_MinMaxs( property_number).GetMin());
          // write the name of the property out (every 4th value in the command line flag)
          m_TableOutputFile
              << '\n' << m_AveMinMaxStdPropertyFlag->GetParameterList()( property_number)->GetValue()
              << '\t' << m_MeanStds( property_number).GetAverage()
              << '\t' << stdev
              << '\t' << m_MinMaxs( property_number).GetMin()
              << '\t' << m_MinMaxs( property_number).GetMax()
              << '\t' << range
              << '\t' << ( stdev > 0.0 ? range / stdev : 0.0);
        }
        io::File::CloseClearFStream( m_TableOutputFile);
      }

      // end
      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief format a string for writing into a csv file
    //! @param STRING string to format
    //! @return the formatted string
    std::string MoleculeProperties::FormatCSVString( const std::string &STRING)
    {
      // if there are no commas, double quotes or new lines, return the string
      if( STRING.find_first_of( ",\"\n") == std::string::npos)
      {
        return STRING;
      }

      // quote the string
      std::string new_string( 1, '"');
      new_string += STRING;
      new_string += '"';

      // if there are no double quotes, return the quoted string
      if( STRING.find( '"') == std::string::npos)
      {
        return new_string;
      }

      // replace " with """
      util::StringReplacement quotes_escaper( util::StringReplacement::e_Any, "\"", std::string( 3, '"'));
      quotes_escaper.ReplaceAllIn( new_string);
      return new_string;
    }

    //! @brief initialize the histograms and properties from the command line flags
    void MoleculeProperties::Initialize() const
    {
      // reset the variables used to hold numeric properties and histograms
      m_NumericProperties.Reset();
      m_NumericHistograms.Reset();
      m_StringProperties.Reset();
      m_TableProperties.Reset();
      m_PropertiesToRemove.Reset();
      m_PropertiesToRename.Reset();
      m_PropertiesToAdd.Reset();
      m_MeanStds.Reset();
      m_MinMaxs.Reset();
      m_StatisticsProperties.Reset();
      m_AddIndexProperty = false;

      // recreate the string histogram vector to hold all the string histograms desired
      m_StringHistograms =
        storage::Vector< storage::Map< std::string, size_t> >( m_HistogramStringPropertyFlag->GetSize());

      // store the property strings for the string histograms
      storage::Vector< std::string> string_props( m_HistogramStringPropertyFlag->GetStringList());
      m_StringProperties.AllocateMemory( string_props.GetSize());

      // construct the string properties
      for
      (
        storage::Vector< std::string>::const_iterator itr( string_props.Begin()), itr_end( string_props.End());
        itr != itr_end;
        ++itr
      )
      {
        m_StringProperties.PushBack( StringProperty( *itr));
      }

      storage::Vector< std::string> statistics_props( m_AveMinMaxStdPropertyFlag->GetStringList());
      m_StatisticsProperties.AllocateMemory( statistics_props.GetSize());
      m_MinMaxs.Resize( statistics_props.GetSize());
      m_MeanStds.Resize( statistics_props.GetSize());

      for
      (
        storage::Vector< std::string>::const_iterator itr( statistics_props.Begin()), itr_end( statistics_props.End());
        itr != itr_end;
        ++itr
      )
      {
        // create a label
        const util::ObjectDataLabel label( *itr);

        // if this looks like an atom property, try to make it one
        m_StatisticsProperties.PushBack( descriptor::CheminfoProperty( label));
      }

      // store the property strings for tabular output
      storage::Vector< std::string> table_props( m_TablePropertiesFlag->GetStringList());
      if( !table_props.IsEmpty())
      {
        m_TableProperties.AllocateMemory( table_props.GetSize());

        // open the file
        io::File::MustOpenOFStream( m_TableOutputFile, m_TableOutputFilenameFlag->GetFirstParameter()->GetValue());

        // write out the csv header
        m_TableOutputFile << "Index";

        // construct the string properties
        for
        (
          storage::Vector< std::string>::const_iterator itr( table_props.Begin()), itr_end( table_props.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !util::StartsWith( *itr, "Define("))
          {
            m_TableProperties.PushBack( StringProperty( *itr));
            m_TableOutputFile << ',' << FormatCSVString( *itr);
          }
        }
        m_TableOutputFile << '\n';
      }

      // make properties to access the numeric data used for each numeric property
      for
      (
        size_t property_number( 0), number_properties( m_HistogramNumericPropertyFlag->GetSize());
        property_number < number_properties;
        property_number += 4
      )
      {
        // if the property is already stored on the molecule, then use the stored value
        // otherwise, it will need to be calculated (if possible), so make a property to calculate
        // the values of this property in case it has not already been stored on the small molecules
        const util::ObjectDataLabel data_label
        (
          m_HistogramNumericPropertyFlag->GetParameterList()( property_number)->GetValue()
        );

        m_NumericProperties.PushBack( descriptor::CheminfoProperty( data_label));

        m_NumericHistograms.PushBack
        (
          math::Histogram
          (
            m_HistogramNumericPropertyFlag->GetParameterList()( property_number + 1)->GetNumericalValue< double>(),
            m_HistogramNumericPropertyFlag->GetParameterList()( property_number + 2)->GetNumericalValue< double>(),
            m_HistogramNumericPropertyFlag->GetParameterList()( property_number + 3)->GetNumericalValue< size_t>()
          )
        );
      }

      // setup property renaming
      const storage::Vector< std::string> renamings( m_RenamePropertiesFlag->GetStringList());
      for
      (
        storage::Vector< std::string>::const_iterator
          itr_prop( renamings.Begin()), itr_prop_end( renamings.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        const std::string &current_property_name( *itr_prop);
        ++itr_prop;
        const std::string &new_property_name( *itr_prop);
        m_PropertiesToRename[ current_property_name] = new_property_name;
      }

      // setup property removal
      m_PropertiesToRemove = m_RemovePropertiesFlag->GetStringList();

      // setup property addition
      const storage::Vector< std::string> property_strings( m_AddPropertiesFlag->GetStringList());
      for
      (
        storage::Vector< std::string>::const_iterator
          itr_prop( property_strings.Begin()), itr_prop_end( property_strings.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        if( *itr_prop == "Index") // add an index to each molecule in the ensemble
        {
          m_AddIndexProperty = true;
          continue;
        }

        // create a label for this property
        util::ObjectDataLabel label( *itr_prop);
        std::string prop_name( label.ToString());
        std::stringstream err_stream;

        // If the label is not a string property, then it must be a numeric (atom or small molecule) property
        if( !util::Implementation< chemistry::StringPropertyInterface>::CreateIfPossible( label).second)
        {
          // so make a numeric string property to adapts the atom or small molecule property
          label = util::ObjectDataLabel( "", "Numeric", storage::Vector< util::ObjectDataLabel>( 1, label));
        }

        // create an implementation to hold the function that converts a small molecule to a string
        StringProperty descriptor( label);
        BCL_Assert( descriptor.IsDefined(), *itr_prop + " is not a known property");

        m_PropertiesToAdd.PushBack( storage::Pair< std::string, StringProperty>( prop_name, descriptor));
      }
      if( m_OutputFilenameFlag->GetFlag())
      {
        io::File::MustOpenOFStream( m_MoleculeOutputFile, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      }
    }

    //! @brief add the data from a small molecule to the histograms
    //! @param MOLECULE the small molecule whose data to add
    //! @param MOLECULE_INDEX the index of the molecule currently operated on
    void MoleculeProperties::AddData( chemistry::FragmentComplete &MOLECULE, const size_t &MOLECULE_INDEX) const
    {
      // add properties
      if( m_AddIndexProperty)
      {
        MOLECULE.StoreProperty( "Index", util::Format()( MOLECULE_INDEX));
      }

       // add properties
      for
      (
        storage::List< storage::Pair< std::string, StringProperty> >::const_iterator
          itr_prop( m_PropertiesToAdd.Begin()), itr_prop_end( m_PropertiesToAdd.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        MOLECULE.StoreProperty( itr_prop->First(), itr_prop->Second()->operator()( MOLECULE));
      }

      // add strings as properties
      if( m_AddPropertyStringsFlag->GetFlag())
      {
        for
        (
            size_t property_number( 0), number_properties( m_AddPropertyStringsFlag->GetSize());
            property_number < number_properties - 1;
            property_number += 2
        )
        {
          const std::string label( m_AddPropertyStringsFlag->GetParameterList()( property_number)->GetValue());
          const std::string value( m_AddPropertyStringsFlag->GetParameterList()( property_number + 1)->GetValue());
          MOLECULE.StoreProperty( label, value);
        }
      }

      // rename properties
      for
      (
        storage::Map< std::string, std::string>::const_iterator
          itr_prop( m_PropertiesToRename.Begin()), itr_prop_end( m_PropertiesToRename.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        const std::string &current_property_name( itr_prop->first);
        if( MOLECULE.IsPropertyStored( current_property_name))
        {
          MOLECULE.StoreProperty( itr_prop->second, MOLECULE.GetMDLProperty( current_property_name));
          MOLECULE.RemoveProperty( current_property_name);
        }
      }

      // add the small molecule's data to all string histograms
      for
      (
        size_t property_number( 0), number_properties( m_StringHistograms.GetSize());
        property_number < number_properties;
        ++property_number
      )
      {
        // get the property using the property retrieval method
        std::string property_value( m_StringProperties( property_number)->operator()( MOLECULE));

        // split the strings retrieved from the property
        storage::Vector< std::string> strings( util::SplitString( property_value));

        // get a reference to the map/histogram
        storage::Map< std::string, size_t> &histogram( m_StringHistograms( property_number));

        // add the strings to the histogram
        for
        (
          storage::Vector< std::string>::const_iterator itr_string( strings.Begin()), itr_string_end( strings.End());
          itr_string != itr_string_end;
          ++itr_string
        )
        {
          if( !itr_string->empty())
          {
            histogram[ *itr_string]++;
          }
        }
      }

      // add data from the small molecule to all numeric histograms
      for
      (
        size_t property_number( 0), number_properties( m_NumericHistograms.GetSize());
        property_number < number_properties;
        ++property_number
      )
      {
        // get a reference to the actual property
        descriptor::Base< chemistry::AtomConformationalInterface, float> &property
        (
          *m_NumericProperties( property_number)
        );

        // get a reference to the histogram
        math::Histogram &numeric_histogram( m_NumericHistograms( property_number));

        property.SetObject( MOLECULE);

        // walk over the iterator
        for
        (
          descriptor::Iterator< chemistry::AtomConformationalInterface> itr( property.GetType(), MOLECULE);
          itr.NotAtEnd();
          ++itr
        )
        {
          // calculate the values of the property, if possible
          linal::VectorConstReference< float> values( property( itr));

          // compute a histogram of the values found for this molecule, then append the results
          // to the property_histogram
          for
          (
            linal::VectorConstReference< float>::const_iterator
              itr_values( values.Begin()), itr_values_end( values.End());
            itr_values != itr_values_end;
            ++itr_values
          )
          {
            numeric_histogram.PushBack( *itr_values);
          }
        }
      }

      // add data to the statistics
      size_t statistic_property_number( 0);
      for
      (
        storage::Vector< descriptor::CheminfoProperty>::iterator
          itr( m_StatisticsProperties.Begin()), itr_end( m_StatisticsProperties.End());
        itr != itr_end;
        ++itr, ++statistic_property_number
      )
      {
        // get a reference to the actual property
        descriptor::Base< chemistry::AtomConformationalInterface, float> &property( **itr);
        property.SetObject( MOLECULE);

        math::RunningAverageSD< float> &ave_std( m_MeanStds( statistic_property_number));
        math::RunningMinMax< float> &min_max( m_MinMaxs( statistic_property_number));

        // walk over the iterator
        for
        (
          descriptor::Iterator< chemistry::AtomConformationalInterface> itr_desc( property.GetType(), MOLECULE);
          itr_desc.NotAtEnd();
          ++itr_desc
        )
        {
          // calculate the values of the property, if possible
          linal::VectorConstReference< float> values( property( itr_desc));

          // compute a histogram of the values found for this molecule, then append the results
          // to the property_histogram
          for
          (
            linal::VectorConstReference< float>::const_iterator
              itr_values( values.Begin()), itr_values_end( values.End());
            itr_values != itr_values_end;
            ++itr_values
          )
          {
            ave_std += *itr_values;
            min_max += *itr_values;
          }
        }
      }

      if( !m_TableProperties.IsEmpty())
      {
        // write out the index of the molecule
        m_TableOutputFile << MOLECULE_INDEX;

        // write out csv data to the table
        for
        (
          storage::Vector< StringProperty>::const_iterator
            itr( m_TableProperties.Begin()), itr_end( m_TableProperties.End());
          itr != itr_end;
          ++itr
        )
        {
          // write the property out to the file
          m_TableOutputFile << ',' << FormatCSVString( ( **itr)( MOLECULE));
        }
        m_TableOutputFile << '\n';
      }

      if( m_OutputFilenameFlag->GetFlag())
      {
        if( m_RemoveAllPropertiesFlag->GetFlag())
        {
          MOLECULE.ResetStoredProperties();
          MOLECULE.WriteMDL( m_MoleculeOutputFile);
        }
        else
        {
          for
          (
              storage::Vector< std::string>::const_iterator
              itr_prop( m_PropertiesToRemove.Begin()), itr_prop_end( m_PropertiesToRemove.End());
              itr_prop != itr_prop_end;
              ++itr_prop
          )
          {
            MOLECULE.RemoveProperty( *itr_prop);
          }
          MOLECULE.WriteMDL( m_MoleculeOutputFile);
        }
      }

    }

    const ApplicationType MoleculeProperties::MoleculeProperties_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeProperties(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl

