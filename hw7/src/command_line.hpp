#ifndef PARSER_HPP
#define PARSER_HPP

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <stack>
#include <queue>
#include <algorithm>

#include "commands.hpp"

using namespace std;

static const char quicksave_dir[] = "data/quicksave.scn";

struct Line {
private:
    Line(); // default constructor disabled

public:
    static const int cmd_line_buffer_len = 256;

    char untokenized_argv[cmd_line_buffer_len];
    char argv[cmd_line_buffer_len];
    vector<char*> tokens;

    Line(char* argv);
    ~Line();

    const int toCommandID() const;
    const Command& toCommand() const;
};

class CommandLine {
private:

    static CommandLine* cmd_line;

    bool running;
    stack<Line*> state;

    vector<const Line*> history;

    explicit CommandLine();
    ~CommandLine();

public:
    static void init();
    static bool active();
    static void pause();
    static void run();
    static void kill();

    static const Line* getState();
    static void clearState();

    static const vector<const Line*>& getHistory();
    static void clearHistory();

    static void readLine(istream& input);
};

#endif