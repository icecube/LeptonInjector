#include <iostream>
#include <stdexcept>
#include <unistd.h>
#include <map>

#include "inc/tools.h"

// ---------------TEST FUNCTION BORROWED FROM PHOTOSPLINE-----------------------------------
// https://github.com/IceCubeOpenSource/photospline
// Published under the BSD 2-Clause "Simplified License" with notice below
/*
Copyright (c) 2009-2012, Nathan Whitehorn and Jakob van Santen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


struct test_exception : public std::runtime_error{
	test_exception(const std::string& msg):std::runtime_error(msg){}
};

void emit_error(const std::string& file, size_t line,
				const std::string& criterion, const std::string& message){
	std::ostringstream ss;
	ss << file << ':' << line << "\n\t";
	if(message.empty())
		ss << "Assertion failed: \n";
	else
		ss << message << ": \n";
	ss << '\t' << criterion << std::endl;
	throw test_exception(ss.str());
}

std::map<std::string,void(*)()>&
test_registry()
{
	static std::map<std::string,void(*)()> *registry = new std::map<std::string,void(*)()>;
	return *registry;
}

int main(int argc, char* argv[]){
	for(int i=1; i<argc; i++){
		std::string arg=argv[i];
		if(arg=="WORKING_DIRECTORY"){
			if(i+1>=argc){
				std::cerr << "WORKING_DIRECTORY not specified" << std::endl;
				return(1);
			}
			if(chdir(argv[i+1])!=0){
				std::cerr << "Failed to change working directory to " << argv[i+1] << std::endl;
				return(1);
			}
			i++;
		}
	}
	
	std::cout << "Running " << test_registry().size() << " tests" << std::endl;
	bool all_pass=true;
	size_t passes=0, failures=0;
	for(std::map<std::string,void(*)()>::const_iterator test=test_registry().begin();
		test!=test_registry().end(); test++){
		bool pass=false;
		std::cout << test->first << ": ";
		std::cout.flush();
		try{
			(test->second)();
			pass=true;
		}catch(test_exception& ex){
			std::cout << "FAIL\n " << ex.what() << std::endl;
		}catch(std::exception& ex){
			std::cout << "FAIL\n Exception: " << ex.what() << std::endl;
		}catch(...){
			std::cout << "FAIL\n Unknown object thrown" << std::endl;
		}
		if(pass)
			std::cout << "PASS" << std::endl;
		(pass?passes:failures)++;
		all_pass &= pass;
	}
	std::cout << passes << " test" << (passes!=1?"s":"") << " pass"
	<< (passes!=1?"":"es") << ", "
	<< failures << " fail" << (failures!=1?"":"s") << std::endl;
	return(all_pass ? 0 : 1);
}