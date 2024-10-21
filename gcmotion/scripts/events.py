r"""
Events are used by SciPy's ``solve_ivp`` solver to locate where certain conditions
are met during the solving, and optionally terminate the solving after a certain
amount of times those conditions are found. More info 
`here <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_.

Events are *function handles* of **analytical** funtions, and they trigger when their
``return`` expression is equal to 0. Since ``solve_ivp`` uses a root find algorithm to 
locate the exact position of the event, the expressions must be 0 **and** change sings 
around the root.

We can use multiple events, and the solver will return a list with the triggered positions
for each of them.

Example
-------

This event triggers whenever :math:`\theta` takes the value "theta0" and terminates 
after "terminal" triggers, which we can pass as arguments. The ``direction``
parameter only checks the **sign** of the direction.

.. code-block:: python

    def when_theta(theta0, terminal):

        def event(t, S):
            return S[0] - theta0

        event.terminal = terminal
        event.direction = 1

        return event

.. rubric:: Availiable events
    :heading-level: 4

"""


def when_theta(theta0, terminal):
    r"""
    Triggers when :math:`\theta` is equal to ``theta0`` and
    terminates after ``terminal`` times (Starting position included).

    Setting ``terminal=0`` makes the event non-terminal.

    Parameters
    ----------
    theta0 : int | float
        The :math:`\theta` value that triggers the event.
    terminal : int
        The event's termination number.

    Returns
    -------

    event : function handle
        The event function handle to be passed to solve_ivp.
    """

    def event(t, S):
        return S[0] - theta0

    if terminal != 0:
        event.terminal = terminal

    return event


def when_psi(psi0, terminal):
    r"""
    Triggers when :math:`\psi` is equal to ``psi0`` and
    terminates after ``terminal`` times (Starting position included).

    .. note::
        This event is used to stop the solver from calculating
        unecessarily long orbits, when 1 should be enough. Since :math:`\psi`
        seems to mostly be a well behaved periodic function and stays bounded
        even for the nastiest initial condition, its a good enough critirion
        for stopping the solver. In case it fails, the solver will just keep
        going, so no harm is done.

    Setting ``terminal=0`` makes the event non-terminal.

    Parameters
    ----------
    theta0 : int | float
        The :math:`\theta` value that triggers the event.
    terminal : int
        The event's termination number.

    Returns
    -------

    event : function handle
        The event function handle to be passed to solve_ivp.
    """

    def event(t, S):
        return S[1] - psi0

    if terminal != 0:
        event.terminal = terminal

    return event
