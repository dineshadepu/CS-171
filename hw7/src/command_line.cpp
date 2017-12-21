#include "command_line.hpp"

/*********************************** Line *************************************/

Line::Line(char* argv) : untokenized_argv(), argv(), tokens() {
    assert(strlen(argv) <= cmd_line_buffer_len);

    // set argv
    strcpy(this->untokenized_argv, argv);
    strcpy(this->argv, argv);
    
    // set tokens
    char* tok = strtok(this->argv, " ");
    assert(tok != NULL);
    while (tok) {
        this->tokens.push_back(tok);
        tok = strtok(NULL, " ");
    }
}
Line::~Line() {

}

const int Line::toCommandID() const {
    // take first arg as the command name
    CommandName cmd_name(this->tokens[0]);

    unordered_map<CommandName, int, CommandNameHasher>::const_iterator
        cmd_id_it = Commands::name_to_id.find(cmd_name);

    // if command can't be found return -1 (invalid_cmd_id) else return found id
    if (cmd_id_it == Commands::name_to_id.end()) {
        return -1;
    } else {
        return cmd_id_it->second;
    }
}

const Command& Line::toCommand() const {
    // take first arg as the command name
    CommandName cmd_name(this->tokens[0]);

    unordered_map<CommandName, int, CommandNameHasher>::const_iterator
        cmd_id_it = Commands::name_to_id.find(cmd_name);

    // if command can't be found print error message and return invalid_cmd
    if (cmd_id_it == Commands::name_to_id.end()) {
        fprintf(stderr, "ERROR invalid command %s\n", cmd_name.name);
        return Commands::id_to_cmd[-1];
    } else {
        // get valid command
        unordered_map<int, Command>::const_iterator cmd_it =
            Commands::id_to_cmd.find(cmd_id_it->second);
        assert(cmd_it != Commands::id_to_cmd.end());

        // if expected_argc range does not match the given argc, print error
        // message and return invalid_cmd
        if (cmd_it->second.min_expected_argc > this->tokens.size() ||
            cmd_it->second.max_expected_argc < this->tokens.size())
        {
            fprintf(stderr, "ERROR invalid argc for command %s. expected values in the range [%d, %d] received %d\n",
                cmd_name.name,
                cmd_it->second.min_expected_argc,
                cmd_it->second.max_expected_argc,
                (int) this->tokens.size());
            return Commands::id_to_cmd[-1];
        } else { // else Line is valid so return found Command
            return cmd_it->second;
        }
    }
}

/******************************** CommandLine *********************************/

CommandLine* CommandLine::cmd_line = NULL;

CommandLine::CommandLine() : running(true), state(), history() {

}
CommandLine::~CommandLine() {

}

void CommandLine::init() {
    Commands::init();
    if (!cmd_line) {
        cmd_line = new CommandLine();
        printf("CS 171 Modeling Language\n");
        printf("type help for the instruction manual\n");
    } else {
        printf("WARNING CommandLine::cmd_line already exists. skipping.\n");
    }
}

void CommandLine::kill() {
    delete cmd_line;
    exit(0);
}

bool CommandLine::active() {
    if (cmd_line) {
        return cmd_line->running;
    }
    fprintf(stderr, "ERROR CommandLine::instance CommandLine has not been initialized yet\n");
    exit(1);
}
void CommandLine::pause() {
    if (cmd_line) {
        cmd_line->running = false;
    } else {
        fprintf(stderr, "ERROR CommandLine::instance CommandLine has not been initialized yet\n");
        exit(1);
    }
}
void CommandLine::run() {
    if (cmd_line) {
        cmd_line->running = true;
    } else {
        fprintf(stderr, "ERROR CommandLine::instance CommandLine has not been initialized yet\n");
        exit(1);
    }
}

const Line* CommandLine::getState() {
    if (cmd_line) {
        if (cmd_line->state.size() > 0) {
            return cmd_line->state.top();
        }
        return NULL;
    }
    fprintf(stderr, "ERROR CommandLine::instance CommandLine has not been initialized yet\n");
    exit(1);
}

void CommandLine::clearState() {
    if (cmd_line->state.size() > 0) {
        cmd_line->state.pop();
    }
}

const vector<const Line*>& CommandLine::getHistory() {
    if (cmd_line) {
        return cmd_line->history;
    }
    fprintf(stderr, "ERROR CommandLine::instance CommandLine has not been initialized yet\n");
    exit(1);
}

void CommandLine::clearHistory() {
    if (cmd_line) {
        cmd_line->history.clear();
    }
}

void CommandLine::readLine(istream& input) {
    assert(cmd_line);

    // read in new line
    char new_line_buffer[Line::cmd_line_buffer_len];
    input.getline(new_line_buffer, Line::cmd_line_buffer_len);

    // find line start based on where the first non space char is
    int line_begin = 0;
    char* blank_char = strchr(new_line_buffer, ' ');
    while (blank_char != NULL && line_begin == blank_char - new_line_buffer) {
        line_begin++;
        blank_char = strchr(blank_char + 1, ' ');
    }

    // find line end after removing any comments
    int line_end = strlen(new_line_buffer);
    char* comment_start = strchr(new_line_buffer, '#');
    if (comment_start) {
        *comment_start = '\0';
        line_end = comment_start - new_line_buffer;
    }

    if (line_begin < line_end) {
        char trimmed_new_line_buffer[Line::cmd_line_buffer_len];
        memcpy(
            trimmed_new_line_buffer,
            new_line_buffer + line_begin,
            line_end - line_begin + 1);

        // setup a new Line instance
        Line* new_line = new Line(trimmed_new_line_buffer);

        // if valid execute (invalid commands will return invalid_cmd which has
        // a NULL pointer for action)
        const Command& cmd = new_line->toCommand();
        if (cmd.action) {
            // if cmd is a state impacting command, update state
            if (cmd.type == STATE) {
                clearState();
                cmd_line->state.push(new_line);
            }

            bool success =
                cmd.action(new_line->tokens.size(), new_line->tokens.data());

            if (success) {
                cmd_line->history.push_back(new_line);
            } else {
                delete new_line;
            }
        }
    }
}

