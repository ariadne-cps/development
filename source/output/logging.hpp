/***************************************************************************
 *            output/logging.hpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file output/logging.hpp
 *  \brief Support for writing debugging output to a logging stream.
 */

#ifndef ARIADNE_LOGGING_HPP
#define ARIADNE_LOGGING_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <utility>
#include <vector>
#include <queue>
#include <map>
#include <sstream>
#include "utility/macros.hpp"
#include "utility/pointer.hpp"

// Set the verbosity, i.e., the maximum level to be shown from println calls
#define ARIADNE_LOG_SET_VERBOSITY(level) Logger::instance().configuration().set_verbosity(level);
// Automatic level increase/decrease in a scope; meant to be used once within a function, at top scope; necessary for print holding.
#define ARIADNE_LOG_SCOPE_CREATE auto logscopemanager = LogScopeManager(ARIADNE_PRETTY_FUNCTION);
// Managed level increase/decrease around the function fn; if the function throws, manual decrease of the proper level is required.
#define ARIADNE_LOG_RUN_AT(level,fn) Logger::instance().increase_level(level); fn; Logger::instance().decrease_level(level);
// Mute the logger for the function fn; if the function throws, manual decrease of the proper level is required.
#define ARIADNE_LOG_RUN_MUTED(fn) Logger::instance().mute_increase_level(); fn; Logger::instance().mute_decrease_level();
// Print one line at the current level; the text shouldn't have carriage returns, but for efficiency purposes this is not checked.
#define ARIADNE_LOG_PRINTLN(text) { if (!Logger::instance().is_muted_at(0)) { std::ostringstream logger_stream; logger_stream << std::boolalpha << text; Logger::instance().println(0,logger_stream.str()); } }
// Print one line at an increased level with respect to the current one; the text shouldn't have carriage returns, but for efficiency purposes this is not checked.
#define ARIADNE_LOG_PRINTLN_AT(level,text) { if (!Logger::instance().is_muted_at(level)) { std::ostringstream logger_stream; logger_stream << std::boolalpha << text; Logger::instance().println(level,logger_stream.str()); } }
// Print variable in one line at the current level, using the formatting convention.
#define ARIADNE_LOG_PRINTLN_VAR(var) { if (!Logger::instance().is_muted_at(0)) { std::ostringstream logger_stream; logger_stream << std::boolalpha << #var << " = " << var; Logger::instance().println(0,logger_stream.str()); } }
// Print variable in one line at the increased level with respect to the current one, using the formatting convention.
#define ARIADNE_LOG_PRINTLN_VAR_AT(level,var) { if (!Logger::instance().is_muted_at(level)) { std::ostringstream logger_stream; logger_stream << std::boolalpha << #var << " = " << var; Logger::instance().println(level,logger_stream.str()); } }
// Print a text at the bottom line, holding it until the function scope ends; this requires creation of the scope.
// Nested calls in separate functions append to the held line.
// The text for obvious reasons shouldn't have newlines and carriage returns; for efficiency purposes this is not checked.
#define ARIADNE_LOG_SCOPE_PRINTHOLD(text) { if (!Logger::instance().is_muted_at(0)) { std::ostringstream logger_stream; logger_stream << std::boolalpha << text; Logger::instance().hold(logscopemanager.scope(),logger_stream.str()); } }

namespace Ariadne {

class LoggableSmartThread;

inline unsigned int get_verbosity(int argc, const char* argv[]) {
    if(argc>1) {
        if(std::strcmp(argv[1],"-v")==0) {
            if(argc>2) {
                int val = std::atoi(argv[2]);
                if (val < 0) std::cerr << "Verbosity should be a non-negative value.\n";
                return static_cast<unsigned int>(val);
            }
        } else {
            std::cerr << "Unrecognised command-line option \"" << argv[1] << "\"\n";
        }
    }
    return 0;
}

//! \brief A styling for a text character
//! \details Refer to https://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html#256-colors
//! for the detailed 256 map of colors for text and background
struct TerminalTextStyle {
    TerminalTextStyle(uint8_t fontcolor, uint8_t bgcolor, bool bold, bool underline);

    //! \brief Return the string code that applies the style
    std::string operator()() const;

    //! \brief If a style is actually set
    bool is_styled() const;

    uint8_t fontcolor; // The color of the character
    uint8_t bgcolor; // The color of the background of the character box
    bool bold; // If the character must be highlighted in bold
    bool underline; // If the character must be underlined

    inline static std::string RESET = "\u001b[0m";
};

static TerminalTextStyle TT_STYLE_NONE(0,0,false,false);
static TerminalTextStyle TT_STYLE_DARKBROWN(95,0,false,false);
static TerminalTextStyle TT_STYLE_DARKORANGE(130,0,false,false);
static TerminalTextStyle TT_STYLE_LIGHTBROWN(136,0,false,false);
static TerminalTextStyle TT_STYLE_BRIGHTORANGE(208,0,false,false);
static TerminalTextStyle TT_STYLE_CREAM(223,0,false,false);
static TerminalTextStyle TT_STYLE_OBSIDIAN(237,0,false,false);
static TerminalTextStyle TT_STYLE_DARKGREY(244,0,false,false);
static TerminalTextStyle TT_STYLE_LIGHTGREY(251,0,false,false);

//! \brief A set of text styles for different elements of the logging interface
struct TerminalTextTheme {
    TerminalTextTheme();
    TerminalTextTheme(TerminalTextStyle level_number, TerminalTextStyle level_shown_separator, TerminalTextStyle level_hidden_separator,
                      TerminalTextStyle multiline_separator, TerminalTextStyle assignment_comparison, TerminalTextStyle miscellaneous_operator,
                      TerminalTextStyle round_parentheses, TerminalTextStyle square_parentheses, TerminalTextStyle curly_parentheses,
                      TerminalTextStyle colon, TerminalTextStyle comma, TerminalTextStyle number, TerminalTextStyle at, TerminalTextStyle keyword);

    TerminalTextStyle level_number; // The number representing the log level
    TerminalTextStyle level_shown_separator; // The separator when the level is shown
    TerminalTextStyle level_hidden_separator; // The separator when the level is hidden

    TerminalTextStyle multiline_separator; // The separator for multiple lines output
    TerminalTextStyle assignment_comparison; // Characters: =, >, <, !
    TerminalTextStyle miscellaneous_operator; // Characters: +, -, *, /, \, ^, |, &, %

    TerminalTextStyle round_parentheses; // Characters: ( )
    TerminalTextStyle square_parentheses; // Characters: [ ]
    TerminalTextStyle curly_parentheses; // Characters: { }

    TerminalTextStyle colon; // Character: :
    TerminalTextStyle comma; // Character: ,
    TerminalTextStyle number; // A 0-9 digit, or a dot after a digit, but not if parsing a numbered_variable

    TerminalTextStyle at; // Character: @
    TerminalTextStyle keyword; // Specific words, not pre/post-fixed by an alphanumeric character: virtual, const, true, false, inf

    bool has_style() const; // Whether at least one style is set

    friend OutputStream& operator<<(OutputStream& os, TerminalTextTheme const& theme);
};

//! \brief Empty theme, for not forcing any theme
static TerminalTextTheme TT_THEME_NONE = TerminalTextTheme();
//! \brief Theme for black background
static TerminalTextTheme TT_THEME_DARK = TerminalTextTheme(TT_STYLE_BRIGHTORANGE,TT_STYLE_BRIGHTORANGE,TT_STYLE_DARKGREY,
                                                           TT_STYLE_DARKGREY,TT_STYLE_BRIGHTORANGE,TT_STYLE_LIGHTBROWN,
                                                           TT_STYLE_DARKORANGE,TT_STYLE_DARKGREY,TT_STYLE_LIGHTBROWN,
                                                           TT_STYLE_LIGHTBROWN,TT_STYLE_DARKGREY,TT_STYLE_CREAM,
                                                           TT_STYLE_DARKGREY,TT_STYLE_DARKGREY);
//! \brief Theme for white background
static TerminalTextTheme TT_THEME_LIGHT = TerminalTextTheme(TT_STYLE_BRIGHTORANGE,TT_STYLE_BRIGHTORANGE,TT_STYLE_DARKGREY,
                                                            TT_STYLE_DARKGREY,TT_STYLE_DARKBROWN,TT_STYLE_DARKBROWN,
                                                            TT_STYLE_DARKGREY,TT_STYLE_OBSIDIAN,TT_STYLE_LIGHTBROWN,
                                                            TT_STYLE_LIGHTBROWN,TT_STYLE_OBSIDIAN,TT_STYLE_DARKORANGE,
                                                            TT_STYLE_DARKGREY,TT_STYLE_DARKGREY);

//! \brief Support class for managing log level increase/decrease in a scope
//! \details Since it is possible to capture the scope string of a function only, nested scopes in a function have the same scope string
//! (but work as expected in terms of level management)
class LogScopeManager {
  public:
    LogScopeManager(std::string scope);
    std::string scope() const;
  public:
    virtual ~LogScopeManager();
  private:
    mutable std::string _scope;
};

//! \brief Supported kinds of log messages
enum RawMessageKind { PRINTLN, HOLD, RELEASE };

//! \brief A log message in raw form, as submitted to a concurrent scheduler
//! \details This does not hold any thread identifier information yet, for efficiency
struct LogThinRawMessage {

    LogThinRawMessage(std::string scope, unsigned int level, std::string text);

    std::string scope;
    unsigned int level;
    std::string text;

    RawMessageKind kind() const;
};

//! \brief Full log message information in raw form, before formatting for actual output
struct LogRawMessage : public LogThinRawMessage {
    LogRawMessage(std::string id, LogThinRawMessage msg) : LogThinRawMessage(msg), identifier(id) { }
    LogRawMessage(std::string id, std::string scp, unsigned int lvl, std::string txt) : LogRawMessage(id,LogThinRawMessage(scp,lvl,txt)) { }
    LogRawMessage(std::string scp, unsigned int lvl, std::string txt) : LogThinRawMessage(scp,lvl,txt), identifier("") { }
    std::string identifier;
};

class LoggerData;

//! \brief The policy for printing the thread id with respect to the log level: NEVER, BEFORE or AFTER
enum class ThreadNamePrintingPolicy { NEVER, BEFORE, AFTER };

OutputStream& operator<<(OutputStream& os, const ThreadNamePrintingPolicy& p);

//! \brief Configuration of visualisation settings for a Logger
class LoggerConfiguration {
  public:

    LoggerConfiguration();

    //! \brief Configuration setters

    //! \brief The depth of printing, where v=0 prevents all printing
    void set_verbosity(unsigned int v);
    //! \brief If true, indents each line with a number of spaces equal to the level
    void set_indents_based_on_level(bool b);
    //! \brief If false, the level is shown on each first line of a print; if true, only when a print has a different level
    void set_prints_level_on_change_only(bool b);
    //! \brief If true, each entrance in a defined scope is logged
    void set_prints_scope_entrance(bool b);
    //! \brief If true, each exit from a defined scope is logged
    void set_prints_scope_exit(bool b);
    //! \brief If true, text that would be printed on multiple lines is correctly preambled on each line with level and indentation
    void set_handles_multiline_output(bool b);
    //! \brief If true, discards from printing all newlines and their immediately following whitespaces
    void set_discards_newlines_and_indentation(bool b);
    //! \brief Decides if and where to append the thread name to the log level (if the latter is shown)
    //! \details The policy is implicitly NEVER if the scheduler is immediate, or only one thread is registered (logging thread excluded)
    void set_thread_name_printing_policy(ThreadNamePrintingPolicy p);

    //! \brief Configuration getters

    unsigned int verbosity() const;
    bool indents_based_on_level() const;
    bool prints_level_on_change_only() const;
    bool prints_scope_entrance() const;
    bool prints_scope_exit() const;
    bool handles_multiline_output() const;
    bool discards_newlines_and_indentation() const;
    ThreadNamePrintingPolicy thread_name_printing_policy() const;

    //! \brief Style theme for terminal output
    void set_theme(TerminalTextTheme const& theme);
    //! \brief Get the current theme used
    TerminalTextTheme const& theme() const;
    //! \brief Add a keyword to the default ones offered in the them, forcing a given style
    //! \details Adding an existing keyword has no effect
    void add_custom_keyword(std::string const& text, TerminalTextStyle const& style);
    //! \brief Add a keyword to the default ones offered in the them, using the default style
    //! \details Adding an existing keyword has no effect
    void add_custom_keyword(std::string const& text);
    //! \brief Get the map of keywords (a TerminalTextStyle equal to TT_STYLE_NONE implies no custom style forced)
    std::map<std::string,TerminalTextStyle> const& custom_keywords() const;

    friend OutputStream& operator<<(OutputStream& os, LoggerConfiguration const& configuration);

  private:
    unsigned int _verbosity;
    bool _indents_based_on_level;
    bool _prints_level_on_change_only;
    bool _prints_scope_entrance;
    bool _prints_scope_exit;
    bool _handles_multiline_output;
    bool _discards_newlines_and_indentation;
    ThreadNamePrintingPolicy _thread_name_printing_policy;

    TerminalTextTheme _theme;
    std::map<std::string,TerminalTextStyle> _custom_keywords;
};

class LoggerSchedulerInterface;

//! \brief A static class for log output handling.
//! Configuration and final printing is done here, while scheduling is
//! done in a separate class.
class Logger {
  private:
    friend class ImmediateLoggerScheduler;
    friend class BlockingLoggerScheduler;
    friend class NonblockingLoggerScheduler;

    Logger();
  public:
    Logger(Logger const&) = delete;
    void operator=(Logger const&) = delete;

    static Logger& instance() {
        static Logger instance;
        return instance;
    }

    void use_immediate_scheduler();
    void use_blocking_scheduler();
    void use_nonblocking_scheduler();

    void register_thread(LoggableSmartThread const& thread);
    void unregister_thread(LoggableSmartThread const& thread);

    void println(unsigned int level_increase, std::string text);
    void hold(std::string scope, std::string text);
    void release(std::string scope);

    void increase_level(unsigned int i);
    void decrease_level(unsigned int i);
    void mute_increase_level();
    void mute_decrease_level();

    bool is_muted_at(unsigned int i) const;

    unsigned int current_level() const;
    std::string current_thread_name() const;
    unsigned int cached_last_printed_level() const;
    std::string cached_last_printed_thread_name() const;

    LoggerConfiguration& configuration();

  private:
    std::string _apply_theme(std::string const& text) const;
    std::string _apply_theme_for_keywords(std::string const& text) const;
    void _print_preamble_for_firstline(unsigned int level, std::string thread_name);
    void _print_preamble_for_extralines(unsigned int level, std::string thread_name);
    std::string _discard_newlines_and_indentation(std::string const& text);
    void _cover_held_columns_with_whitespaces(unsigned int printed_columns);
    void _print_held_line();
    void _println(LogRawMessage const& msg);
    void _hold(LogRawMessage const& msg);
    void _release(LogRawMessage const& msg);
    unsigned int _get_window_columns() const;
    bool _is_holding() const;
    bool _can_print_thread_name() const;
  private:
    static const unsigned int _MUTE_LEVEL_OFFSET = 1024;
    std::vector<LogRawMessage> _current_held_stack;
    unsigned int _cached_num_held_columns;
    unsigned int _cached_last_printed_level;
    std::string _cached_last_printed_thread_name;
    SharedPointer<LoggerSchedulerInterface> _scheduler;
    LoggerConfiguration _configuration;
};

// Global log output file
extern std::ofstream log_file_stream;

//! \brief Redirect logging output to file \a filename.
void redirect_log(const char* filename);

} // namespace Ariadne

#endif // ARIADNE_LOGGING_HPP
