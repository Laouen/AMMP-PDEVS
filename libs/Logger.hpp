/**
 * Copyright (c) 2017, Laouen Mayal Louan Belloli
 * Carleton University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef PMGBP_PDEVS_LOGGER_HPP
#define PMGBP_PDEVS_LOGGER_HPP

// TODO: add option to log to files instead of the standard output
class Logger {
private:
    std::string module_name = "";
public:

    /**
     * @brief Default constructor
     * @details It construct a new instance of Logger without module name.
     */

    Logger() {};

    /**
     * @brief Constructor with module name.
     * @details This constructor takes as parameter a std::string that it is used
     * to set the module name that will be printed with each log.
     *
     * @param other_module_name A std::string with the module name that will
     * use this instance to print logs.
     */
    Logger(std::string other_module_name) {
        module_name = other_module_name;
    }

    /**
     * @brief It set a new module name
     * @details This method allows to set the module name after the instance
     * construction or to modify the current module name of the instance
     *
     * @param other_module_name A std::string with the new module name that will
     * use this instance to print logs
     */
    void setModuleName(std::string other_module_name) {
        module_name = other_module_name;
    }

    /**
     * @brief Prints log messages.
     * @details If #define show_log is not commented, this method prints messages
     * under the [LOG] tag. If #define show_log is commented, thos method does
     * nothing.
     *
     * @param msg The std::string to print
     */
    void log(std::string msg) const {
#ifdef show_log
        std::cout << "[LOG] - ";
		std::cout << "[" + module_name + "] ";
		std::cout << msg << std::endl;
#endif
    }

    /**
     * @brief Prints info messages.
     * @details If #define show_info is not commented, this method prints messages
     * under the [INFO] tag. If #define show_info is commented, thos method does
     * nothing.
     *
     * @param msg The std::string to print
     */
    void info(std::string msg) const {
#ifdef show_info
        std::cout << "[INFO] - ";
		std::cout << "[" + module_name + "] ";
		std::cout << msg << std::endl;
#endif
    }

    /**
     * @brief Prints info messages.
     * @details If #define show_debug is not commented, this method prints messages
     * under the [DEBUG] tag. If #define show_debug is commented, thos method does
     * nothing.
     *
     * @param msg The std::string to print
     */
    void debug(std::string msg) const {
#ifdef show_debug
        std::cout << "[DEBUG] - ";
		std::cout << "[" + module_name + "] ";
		std::cout << msg << std::endl;
#endif
    }

    /**
     * @brief Prints info messages.
     * @details If #define show_error is not commented, this method prints messages
     * under the [ERROR] tag. If #define show_error is commented, thos method does
     * nothing.
     *
     * @param msg The std::string to print
     */
    void error(std::string msg) const {
#ifdef show_error
        std::cout << "[ERROR] - ";
		std::cout << "[" + module_name + "] ";
		std::cout << msg << std::endl;
#endif
    }
};


#endif //PMGBP_PDEVS_LOGGER_HPP
