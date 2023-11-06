% Little code in this part. Just some important tips.

%% Handle vs value classes
% Every class in funkyMage (except fm_eventList) is a "handle" class, which
% means that when you create an object of that class, your variable is a
% reference to that object, but that object is not inherently bound to your
% variable. Assigning the variable to another variable will refer to the
% *same* object. To demonstrate:

% Let's make a task object and fill it in with something
task = fm_task();
task.content = {1, "2 1", "", "", "1 2", "1 2", "1 0"};

% Let's make a copy and empty the content
task2 = task;
task2.content = [];

% Now see if it worked
task2.content % yup, it's empty

% But also
task.content % this is empty too

% This is because when we did 'task2 = task', we merely passed a reference
% to the object, not the object itself.

% Just beware of this. If you want to make a full copy of an object and
% make a new reference to the new object, do:
% task2 = copy(task);
% If you replace 'task2 = task' with this you will see the problem is now
% gone.

% For more, go to:
% https://www.mathworks.com/help/matlab/matlab_oop/comparing-handle-and-value-classes.html

%% Static methods
% If you've looked through some funkyMage help, you'll notice that some
% methods on classes need to be called using a class object
% -- for example:
sim = fm_simulation();
sim.go();

% ...whereas others are just called using the class name e.g.,
funkyMage.addToPath()

% Above, you do not need to create an object from the class 'funkyMage' to
% use the 'addToPath' method. This is because this is a *static* method,
% which means it can be used even without creating an object. Static
% methods can be accessed just through the class name (as seen above).

%% Public vs protected classes
% As you know, classes bundle properties and methods into a template, which
% can be used to create an object. These classes can have several types of
% properties or methods, one important distinction being in their access
% level, which can be 'public', 'protected', or 'private'. If you ever want
% to implement your own features into this toolbox, you can read/write
% existing properties and overridde existing methods as long as they are
% public or protected. If adding features for yourself ever becomes of
% interest, you can look here:
% https://www.mathworks.com/help/simscape/lang/subclassing-and-inheritance.html
% https://www.mathworks.com/help/matlab/matlab_oop/property-attributes.html

%% Contributions
% If you'd like to improve this toolbox, make a contribution through
% GitHub! 