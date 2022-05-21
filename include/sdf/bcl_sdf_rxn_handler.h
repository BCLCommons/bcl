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

#ifndef BCL_SDF_RXN_HANDLER_H_
#define BCL_SDF_RXN_HANDLER_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_mdl_handler.h"
#include "bcl_sdf_molfile_handler.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RXNHandler
    //! @brief takes care of reading/parsing 
    //!        certain line types as well as writing the lines back to create a MDL file from a MDL Handler.
    //!
    //! @see @link example_sdf_rxn_handler.cpp @endlink
    //! @author geanesar, combss, mendenjl
    //! @date Jun 02, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RXNHandler :
      public util::ObjectInterface
    {
    public:

      //! number molecule descriptor lines
      static const size_t s_NumberDescriptionLines;

    private:

    //////////
    // data //
    //////////

      //! rxn description lines
      mutable std::string                     m_Description;

      //! how many reactants are need for this reaction
      mutable size_t                          m_NumberReactants;

      //! how many products are created from this reactions
      mutable size_t                          m_NumberProducts;

      //! The mdl handler for reactants
      mutable storage::Vector< MolfileHandler>   m_ReactantMolfiles;

      //! the mdl handler for products
      mutable storage::Vector< MolfileHandler>   m_ProductMolfiles;

      //! store whether the read lines were a valid reaction
      mutable bool                            m_IsValid;

      //! a mapping of reactive atoms to their atom identity in the reactants
      //! stored as reactive atom (key) to (reactant index,atom number) (value)
      mutable storage::Map< size_t, storage::Pair< size_t, size_t> > m_ReactiveAtomsReactants;

      //! a mapping of reactive atoms to their atom identity in the products
      //! stored as reactive atom (key) to (product index,atom number) (value)
      mutable storage::Map< size_t, storage::Pair< size_t, size_t> > m_ReactiveAtomsProducts;

    public:

      //! s_Instance enables reading/writing ShPtr to this object
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RXNHandler();

      //! @brief constructor with initial input 
      explicit RXNHandler( std::istream &ISTREAM);

      //! @brief Clone function
      //! @return pointer to new RXNHandler
      RXNHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get all description lines stored in handler
      //! @return list of mdl description lines stored in handler
      const std::string &GetDescription() const;

      //! @brief check if the handler is valid
      //! @return true if the mdl handler is in a valid state
      bool IsValid() const;

      //! @brief Get the number of Reactants
      const size_t &GetNumberReactants() const;

      //! @brief Get the number of Reactants
      const size_t &GetNumberProducts() const;

      const storage::Vector< MolfileHandler> &GetReactantHandlers() const
      {
        return m_ReactantMolfiles;
      }

      const storage::Vector< MolfileHandler> &GetProductHandlers() const
      {
        return m_ProductMolfiles;
      }

      //! @brief Get ReactantMdlHandler
      const storage::Vector< MdlHandler> &GetReactantMdlHandlers() const;

      //! @brief Get ProductMdlHandler
      const storage::Vector< MdlHandler> &GetProductMdlHandlers() const;

      //! @brief Get reacting atom mapping on the reactants
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &GetReactiveAtomsInReactants() const;

      //! @brief Get reacting atom mapping on the products
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &GetReactiveAtomsInProducts() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief check if a line is the terminating line
      //! @return true if the line just contains $RXN
      static bool IsRXNDelimiter( const std::string &LINE);

      //! @brief check if a line is the terminating line
      //! @return true if the line just contains $MOL
      static bool IsRXNMolDelimiter( const std::string &LINE);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromRXN( std::istream &ISTREAM);

      //! @brief read RXN from a file buffer as iterators to strings
      //! @param LINE_BEGIN iterator to first line to begin reading (should be "$RXN" but may be blank/whitespace only)
      //! @param LINE_END iterator to one past the last line that should be read
      //! @return iterator to the first line that was not read.  If all were read, should be LINE_END
      storage::List< std::string>::const_iterator ReadFromRXN
      (
        const storage::List< std::string>::const_iterator &LINE_BEGIN,
        const storage::List< std::string>::const_iterator &LINE_END
      );

      static std::ostream &WriteToRXN
      (
        std::ostream &OSTREAM,
        const std::string &DESCRIPTION,
        const storage::Vector< MolfileHandler> &REACTANT_MOLFILES,
        const storage::Vector< MolfileHandler> &PRODUCT_MOLFILES
      );

      //! @brief write a RXN-formatted reaction to an output stream 
      //! @param OSTREAM the output stream
      //! @param DESCRIPTION the description of the reaction
      //! @param REACTANTS the reactants to write
      //! @param PRODUCTS the products to write
      //! @param REACTANT_ATOM_MAP reactive atom map for reactants
      //! @param PRODUCT_ATOM_MAP reactive atom map for products
      //! @param TERMINATION_LINE the line to terminate the RXN with
      //! @return the output stream that was written
      static std::ostream &WriteToRXN
      (
        std::ostream &OSTREAM, 
        const std::string &DESCRIPTION,
        const storage::Vector< chemistry::FragmentComplete> &REACTANTS,
        const storage::Vector< chemistry::FragmentComplete> &PRODUCTS,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTANT_ATOM_MAP,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &PRODUCT_ATOM_MAP
      ); 

      static std::ostream &WriteToRXN
      (
        std::ostream &OSTREAM, 
        const std::string &DESCRIPTION,
        const storage::Vector< storage::Vector< AtomInfo> > &REACTANT_ATOM_INFO,
        const storage::Vector< storage::Vector< BondInfo> > &REACTANT_BOND_INFO,
        const storage::Vector< storage::Vector< AtomInfo> > &PRODUCT_ATOM_INFO,
        const storage::Vector< storage::Vector< BondInfo> > &PRODUCT_BOND_INFO,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTANT_ATOM_MAP,
        const storage::Map< size_t, storage::Pair< size_t, size_t> > &PRODUCT_ATOM_MAP,
        const storage::Vector< std::string> &REACTANT_DESCRIPTIONS = storage::Vector< std::string>(),
        const storage::Vector< std::string> &PRODUCT_DESCRIPTIONS = storage::Vector< std::string>()
      ); 

    protected:

      //! @brief finalize parsing; translates m_Stream into internal members, sets m_Parsed to true, cleans m_Stream
      bool ValidateInput() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set the description; checks # of newlines, ensures that it is exactly s_NumberDescriptionLines
      //! If the # of newlines is < s_NumberDescriptionLines - 1, adds new lines as necessary
      //! Any newlines >= s_NumberDescriptionLines are replaced with spaces
      static std::string StandardizeDescription( const std::string &DESCRIPTION);

      //! @brief test whether a line contains only spaces or is otherwise empty
      //! @param STRING the string to test
      static bool ContainsNonspaceCharacters( const std::string &STRING);

      //! @brief write misc properties as mdl lines into std::ostream
      //! @param OSTREAM ostream to write MiscProperties to
      static void WriteMiscProperties
      (
        std::ostream &OSTREAM,
        const storage::Map< std::string, std::string> &MISC_PROPERTIES
      );

      //! @brief return the data label, if line is a data label line
      //! @param LINE line from mdl section
      //! @return string that constains datalable, string will be empty for non-data lable lines
      static std::string GetMDLDataLabel( const std::string &LINE);

    }; // class RXNHandler

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_RXN_HANDLER_H_

